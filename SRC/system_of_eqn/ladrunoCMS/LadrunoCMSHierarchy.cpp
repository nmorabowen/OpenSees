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

#include "LadrunoCMSHierarchy.h"

#include "LadrunoCMSLocalLanczos.h"
#include "LadrunoCMSMumps.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>

#if defined(_WIN32)
// NOMINMAX: windows.h defines min/max as MACROS, which breaks every std::min /
// std::max in this file (C2589 "illegal token on right side of '::'"). Must be
// defined before windows.h is pulled in.
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")
#elif defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

#include <mpi.h>

#ifdef _WIN32
extern "C" int DSYGVX(
    int *, char *, char *, char *, int *, double *, int *, double *, int *,
    double *, double *, int *, int *, double *, int *, double *, double *, int *,
    double *, int *, int *, int *, int *);
#define LADRUNO_HIERARCHY_DSYGVX DSYGVX
#else
extern "C" int dsygvx_(
    int *, char *, char *, char *, int *, double *, int *, double *, int *,
    double *, double *, int *, int *, double *, int *, double *, double *, int *,
    double *, int *, int *, int *, int *);
#define LADRUNO_HIERARCHY_DSYGVX dsygvx_
#endif

namespace ladruno_cms {
namespace {

struct DenseMatrix {
    int rows = 0;
    int columns = 0;
    std::vector<double> values;

    double &operator()(int row, int column)
    {
        return values[static_cast<std::size_t>(column) * rows + row];
    }
    double operator()(int row, int column) const
    {
        return values[static_cast<std::size_t>(column) * rows + row];
    }
};

enum class KeyKind { Level2Mode, PhysicalEquation, Level1Mode };

struct CoordinateKey {
    KeyKind kind = KeyKind::PhysicalEquation;
    int owner = -1;
    int index = -1;

    bool operator==(const CoordinateKey &other) const
    {
        return kind == other.kind && owner == other.owner && index == other.index;
    }
};

struct CBReduction {
    DenseMatrix transformation;
    SymmetricCSR stiffness;
    SymmetricCSR mass;
    std::vector<int> interior;
    std::vector<int> boundary;
    int retainedModes = 0;
};

struct FineReduction {
    const FineSubdomainInput *input = nullptr;
    CBReduction cb;
    std::vector<CoordinateKey> keys;
};

struct Level1Assembly {
    int coarse = -1;
    std::vector<int> fineIndices;
    std::vector<std::vector<int>> fineToUnique;
    std::vector<CoordinateKey> keys;
    std::vector<int> counts;
    SymmetricCSR stiffness;
    SymmetricCSR mass;
};

struct CoarseReduction {
    int coarse = -1;
    int assemblyIndex = -1;
    CBReduction cb;
    std::vector<CoordinateKey> keys;
};

bool checkedProduct(std::size_t left, std::size_t right, std::size_t &value)
{
    if (left != 0u && right > std::numeric_limits<std::size_t>::max() / left)
        return false;
    value = left * right;
    return true;
}

bool checkedMPIProduct(int left, int right, int &value)
{
    if (left < 0 || right < 0)
        return false;
    const long long product = static_cast<long long>(left) * right;
    if (product > INT_MAX)
        return false;
    value = static_cast<int>(product);
    return true;
}

// ADR-1000 P4 section 1. Peak resident set of THIS process, bytes. The P4 plan
// assumed a Linux CI lane (getrusage ru_maxrss); the harness that actually runs
// CMS today is Windows, so both are wired. Returns 0 when the query fails --
// callers must treat 0 as "unknown", never as "no memory used".
std::size_t peakResidentSetBytes()
{
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS counters;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters)))
        return static_cast<std::size_t>(counters.PeakWorkingSetSize);
    return 0u;
#elif defined(__unix__) || defined(__APPLE__)
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
#if defined(__APPLE__)
        return static_cast<std::size_t>(usage.ru_maxrss);        // bytes
#else
        return static_cast<std::size_t>(usage.ru_maxrss) * 1024u; // kilobytes
#endif
    }
    return 0u;
#else
    return 0u;
#endif
}

DenseMatrix makeDense(int rows, int columns)
{
    DenseMatrix result;
    result.rows = rows;
    result.columns = columns;
    std::size_t size = 0;
    if (rows >= 0 && columns >= 0 &&
        checkedProduct(static_cast<std::size_t>(rows), static_cast<std::size_t>(columns), size))
        result.values.assign(size, 0.0);
    return result;
}

int csrFromDenseUpper(const DenseMatrix &dense, SymmetricCSR &result, std::string &message)
{
    result = SymmetricCSR{};
    std::size_t expected = 0;
    if (dense.rows < 0 || dense.rows != dense.columns ||
        !checkedProduct(static_cast<std::size_t>(dense.rows),
                        static_cast<std::size_t>(dense.columns), expected) ||
        dense.values.size() != expected) {
        message = "cannot convert a non-square dense matrix to symmetric CSR";
        return -1;
    }
    result.dimension = dense.rows;
    result.rowOffsets.push_back(0);
    for (int row = 0; row < dense.rows; ++row) {
        for (int column = row; column < dense.columns; ++column) {
            const double value = dense(row, column);
            if (!std::isfinite(value)) {
                message = "dense hierarchy matrix contains a non-finite value";
                return -2;
            }
            result.columnIndices.push_back(column);
            result.upperValues.push_back(value);
        }
        if (result.columnIndices.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            message = "dense hierarchy matrix exceeds CSR integer range";
            return -3;
        }
        result.rowOffsets.push_back(static_cast<int>(result.columnIndices.size()));
    }
    return result.validate(message);
}

DenseMatrix denseFromCSR(const SymmetricCSR &matrix)
{
    const std::vector<double> rowMajor = matrix.toDense();
    DenseMatrix result = makeDense(matrix.dimension, matrix.dimension);
    if (result.values.empty() && matrix.dimension != 0)
        return DenseMatrix{};
    for (int row = 0; row < matrix.dimension; ++row)
        for (int column = 0; column < matrix.dimension; ++column)
            result(row, column) = rowMajor[static_cast<std::size_t>(row) * matrix.dimension + column];
    return result;
}

DenseMatrix multiply(const DenseMatrix &left, const DenseMatrix &right)
{
    if (left.columns != right.rows)
        return DenseMatrix{};
    DenseMatrix result = makeDense(left.rows, right.columns);
    for (int column = 0; column < right.columns; ++column)
        for (int inner = 0; inner < left.columns; ++inner) {
            const double rightValue = right(inner, column);
            for (int row = 0; row < left.rows; ++row)
                result(row, column) += left(row, inner) * rightValue;
        }
    return result;
}

SymmetricCSR principalSubmatrix(
    const SymmetricCSR &matrix,
    const std::vector<int> &indices,
    std::string &message)
{
    SymmetricCSR result;
    result.dimension = static_cast<int>(indices.size());
    std::vector<int> localOfSource(static_cast<std::size_t>(matrix.dimension), -1);
    for (std::size_t local = 0; local < indices.size(); ++local) {
        const int source = indices[local];
        if (source < 0 || source >= matrix.dimension ||
            localOfSource[static_cast<std::size_t>(source)] >= 0) {
            message = "principal-submatrix indices are invalid or repeated";
            return SymmetricCSR{};
        }
        localOfSource[static_cast<std::size_t>(source)] = static_cast<int>(local);
    }
    std::map<std::pair<int, int>, double> entries;
    for (int sourceRow = 0; sourceRow < matrix.dimension; ++sourceRow) {
        const int localRow = localOfSource[static_cast<std::size_t>(sourceRow)];
        if (localRow < 0)
            continue;
        for (int position = matrix.rowOffsets[static_cast<std::size_t>(sourceRow)];
             position < matrix.rowOffsets[static_cast<std::size_t>(sourceRow + 1)]; ++position) {
            const int sourceColumn = matrix.columnIndices[static_cast<std::size_t>(position)];
            const int localColumn = localOfSource[static_cast<std::size_t>(sourceColumn)];
            if (localColumn < 0)
                continue;
            entries[{std::min(localRow, localColumn), std::max(localRow, localColumn)}] =
                matrix.upperValues[static_cast<std::size_t>(position)];
        }
    }
    result.rowOffsets.push_back(0);
    for (int row = 0; row < result.dimension; ++row) {
        const auto begin = entries.lower_bound({row, row});
        const auto end = entries.upper_bound({row, std::numeric_limits<int>::max()});
        for (auto entry = begin; entry != end; ++entry) {
            result.columnIndices.push_back(entry->first.second);
            result.upperValues.push_back(entry->second);
        }
        result.rowOffsets.push_back(static_cast<int>(result.columnIndices.size()));
    }
    if (result.validate(message) < 0)
        return SymmetricCSR{};
    return result;
}

int crossBlock(
    const SymmetricCSR &matrix,
    const std::vector<int> &rows,
    const std::vector<int> &columns,
    DenseMatrix &result,
    std::string &message)
{
    result = makeDense(static_cast<int>(rows.size()), static_cast<int>(columns.size()));
    std::vector<int> rowMap(static_cast<std::size_t>(matrix.dimension), -1);
    std::vector<int> columnMap(static_cast<std::size_t>(matrix.dimension), -1);
    for (std::size_t local = 0; local < rows.size(); ++local) {
        const int source = rows[local];
        if (source < 0 || source >= matrix.dimension || rowMap[source] >= 0) {
            message = "cross-block row indices are invalid or repeated";
            return -1;
        }
        rowMap[source] = static_cast<int>(local);
    }
    for (std::size_t local = 0; local < columns.size(); ++local) {
        const int source = columns[local];
        if (source < 0 || source >= matrix.dimension || columnMap[source] >= 0) {
            message = "cross-block column indices are invalid or repeated";
            return -2;
        }
        columnMap[source] = static_cast<int>(local);
    }
    for (int sourceRow = 0; sourceRow < matrix.dimension; ++sourceRow) {
        for (int position = matrix.rowOffsets[static_cast<std::size_t>(sourceRow)];
             position < matrix.rowOffsets[static_cast<std::size_t>(sourceRow + 1)]; ++position) {
            const int sourceColumn = matrix.columnIndices[static_cast<std::size_t>(position)];
            const double value = matrix.upperValues[static_cast<std::size_t>(position)];
            if (rowMap[sourceRow] >= 0 && columnMap[sourceColumn] >= 0)
                result(rowMap[sourceRow], columnMap[sourceColumn]) = value;
            if (sourceRow != sourceColumn && rowMap[sourceColumn] >= 0 &&
                columnMap[sourceRow] >= 0)
                result(rowMap[sourceColumn], columnMap[sourceRow]) = value;
        }
    }
    return 0;
}

int congruence(
    const SymmetricCSR &matrix,
    const DenseMatrix &transformation,
    SymmetricCSR &reduced,
    std::string &message)
{
    if (matrix.dimension != transformation.rows) {
        message = "hierarchy congruence dimensions are inconsistent";
        return -1;
    }
    DenseMatrix matrixTimesTransformation =
        makeDense(transformation.rows, transformation.columns);
    for (int column = 0; column < transformation.columns; ++column) {
        std::vector<double> vector(static_cast<std::size_t>(transformation.rows));
        for (int row = 0; row < transformation.rows; ++row)
            vector[static_cast<std::size_t>(row)] = transformation(row, column);
        const std::vector<double> product = matrix.multiply(vector);
        if (product.size() != vector.size()) {
            message = "sparse multiply failed during hierarchy congruence";
            return -2;
        }
        for (int row = 0; row < transformation.rows; ++row)
            matrixTimesTransformation(row, column) = product[static_cast<std::size_t>(row)];
    }
    DenseMatrix denseReduced = makeDense(transformation.columns, transformation.columns);
    for (int column = 0; column < transformation.columns; ++column)
        for (int row = 0; row <= column; ++row) {
            double upper = 0.0;
            double lower = 0.0;
            for (int physical = 0; physical < transformation.rows; ++physical) {
                upper += transformation(physical, row) *
                    matrixTimesTransformation(physical, column);
                lower += transformation(physical, column) *
                    matrixTimesTransformation(physical, row);
            }
            const double value = 0.5 * (upper + lower);
            denseReduced(row, column) = value;
            denseReduced(column, row) = value;
        }
    return csrFromDenseUpper(denseReduced, reduced, message);
}

int reduceCraigBampton(
    const SymmetricCSR &stiffness,
    const SymmetricCSR &mass,
    const std::vector<int> &interior,
    const std::vector<int> &boundary,
    int requestedModes,
    const TwoLevelHierarchyInput &controls,
    int level,
    int subdomain,
    CBReduction &result,
    std::string &message)
{
    result = CBReduction{};
    const int dimension = stiffness.dimension;
    if (mass.dimension != dimension ||
        interior.size() + boundary.size() != static_cast<std::size_t>(dimension)) {
        message = "Craig-Bampton coordinate partition is inconsistent";
        return -1;
    }
    std::vector<bool> seen(static_cast<std::size_t>(dimension), false);
    for (int coordinate : interior) {
        if (coordinate < 0 || coordinate >= dimension || seen[static_cast<std::size_t>(coordinate)]) {
            message = "Craig-Bampton interior contains an invalid coordinate";
            return -2;
        }
        seen[static_cast<std::size_t>(coordinate)] = true;
    }
    for (int coordinate : boundary) {
        if (coordinate < 0 || coordinate >= dimension || seen[static_cast<std::size_t>(coordinate)]) {
            message = "Craig-Bampton boundary contains an invalid coordinate";
            return -3;
        }
        seen[static_cast<std::size_t>(coordinate)] = true;
    }
    result.interior = interior;
    result.boundary = boundary;
    const int numberOfInterior = static_cast<int>(interior.size());
    const int numberOfBoundary = static_cast<int>(boundary.size());
    if (requestedModes < -1) {
        message = "Craig-Bampton mode count must be -1 or nonnegative";
        return -4;
    }

    result.retainedModes = 0;
    DenseMatrix phi = makeDense(numberOfInterior, 0);
    DenseMatrix psi = makeDense(numberOfInterior, numberOfBoundary);
    if (numberOfInterior > 0) {
        const SymmetricCSR stiffnessII = principalSubmatrix(stiffness, interior, message);
        if (stiffnessII.dimension != numberOfInterior)
            return -5;
        MumpsSPD stiffnessFactor;
        if (stiffnessFactor.factorize(stiffnessII, message) < 0) {
            message = "Craig-Bampton K_II is not SPD: " + message;
            return -6;
        }
        if (numberOfBoundary > 0) {
            DenseMatrix stiffnessIB;
            if (crossBlock(stiffness, interior, boundary, stiffnessIB, message) < 0)
                return -7;
            psi.values = stiffnessIB.values;
            for (double &value : psi.values)
                value = -value;
            if (stiffnessFactor.solve(psi.values, numberOfBoundary, message) < 0) {
                message = "Craig-Bampton constraint-mode solve failed: " + message;
                return -7;
            }
        }
        if (requestedModes != 0) {
            const SymmetricCSR massII = principalSubmatrix(mass, interior, message);
            if (massII.dimension != numberOfInterior)
                return -8;
            StaticCondensation condensation;
            if (condenseCoordinateMass(
                    stiffnessII, massII,
                    controls.massRtol, controls.massAtol,
                    condensation, message) < 0) {
                message = "Craig-Bampton interior mass condensation failed: " + message;
                return -9;
            }
            result.retainedModes = requestedModes < 0
                ? condensation.reducedStiffness.dimension
                : std::min(requestedModes, condensation.reducedStiffness.dimension);
            phi = makeDense(numberOfInterior, result.retainedModes);
            LocalLanczosControls lanczos;
            lanczos.tolerance = controls.localTolerance;
            lanczos.maximumOperatorApplications = controls.maximumOperatorApplications;
            lanczos.maximumRestarts = controls.maximumRestarts;
            lanczos.level = level;
            lanczos.subdomain = subdomain;
            LocalEigenResult modes;
            if (solveLocalLanczos(
                    condensation.reducedStiffness,
                    condensation.reducedMass,
                    result.retainedModes,
                    lanczos,
                    modes,
                    message) < 0) {
                message = "Craig-Bampton fixed-interface Lanczos failed: " + message;
                return -11;
            }
            std::vector<double> reconstructed;
            if (reconstructStaticCoordinates(
                    condensation,
                    modes.eigenvectors,
                    result.retainedModes,
                    reconstructed,
                    message) < 0) {
                message = "Craig-Bampton interior mode reconstruction failed: " + message;
                return -12;
            }
            phi.values.swap(reconstructed);
        }
    }

    result.transformation = makeDense(
        dimension, result.retainedModes + numberOfBoundary);
    for (int mode = 0; mode < result.retainedModes; ++mode)
        for (int row = 0; row < numberOfInterior; ++row)
            result.transformation(interior[static_cast<std::size_t>(row)], mode) = phi(row, mode);
    for (int boundaryColumn = 0; boundaryColumn < numberOfBoundary; ++boundaryColumn) {
        const int reducedColumn = result.retainedModes + boundaryColumn;
        for (int row = 0; row < numberOfInterior; ++row)
            result.transformation(interior[static_cast<std::size_t>(row)], reducedColumn) =
                psi(row, boundaryColumn);
        result.transformation(boundary[static_cast<std::size_t>(boundaryColumn)],
                              reducedColumn) = 1.0;
    }
    if (congruence(stiffness, result.transformation, result.stiffness, message) < 0 ||
        congruence(mass, result.transformation, result.mass, message) < 0)
        return -14;
    return 0;
}

struct CoordinateKeyHash {
    std::size_t operator()(const CoordinateKey &key) const noexcept
    {
        std::size_t seed = static_cast<std::size_t>(key.kind);
        seed = seed * 1000003u +
            static_cast<std::size_t>(static_cast<unsigned int>(key.owner));
        seed = seed * 1000003u +
            static_cast<std::size_t>(static_cast<unsigned int>(key.index));
        return seed;
    }
};

// ADR-1000 P4 section 4. This was a linear `std::find` over `unique` for EVERY
// key, i.e. O(totalKeys * uniqueKeys) -- on a Building-1A-sized merge that is
// ~1e9 key comparisons before any arithmetic happens. The hash join is the same
// algorithm with a lookup table; `unique` keeps its first-seen order, which the
// callers depend on (it defines the assembled coordinate numbering).
int compatibilityMaps(
    const std::vector<std::vector<CoordinateKey>> &groups,
    std::vector<std::vector<int>> &maps,
    std::vector<CoordinateKey> &unique,
    std::vector<int> &counts)
{
    maps.clear();
    unique.clear();
    counts.clear();
    std::unordered_map<CoordinateKey, int, CoordinateKeyHash> lookup;
    for (const auto &group : groups) {
        std::vector<int> localMap;
        localMap.reserve(group.size());
        for (const CoordinateKey &key : group) {
            const auto found = lookup.find(key);
            int index = 0;
            if (found == lookup.end()) {
                index = static_cast<int>(unique.size());
                lookup.emplace(key, index);
                unique.push_back(key);
                counts.push_back(0);
            } else {
                index = found->second;
            }
            localMap.push_back(index);
            ++counts[static_cast<std::size_t>(index)];
        }
        maps.push_back(std::move(localMap));
    }
    return unique.empty() ? -1 : 0;
}

int assembleCompatible(
    const std::vector<const SymmetricCSR *> &blocks,
    const std::vector<std::vector<int>> &maps,
    int dimension,
    SymmetricCSR &result,
    std::string &message)
{
    if (blocks.size() != maps.size() || dimension < 1) {
        message = "compatibility assembly dimensions are inconsistent";
        return -1;
    }
    // ADR-1000 P4 section 4 -- the workspace that dominated rank-local memory.
    //
    // This used to materialise a dimension x dimension DENSE buffer, plus a
    // dense copy of every block, and then hand the whole upper triangle to
    // csrFromDenseUpper -- which keeps EVERY entry, zeros included. So the
    // output was structurally dense too: O(dimension^2) paid twice, and a fully
    // dense pattern handed downstream to MUMPS (the same dense pattern that
    // makes MUMPS' auto ordering pick PORD and die in analysis -- see section 21
    // / LEDGER_quirks). At Building-1A scale that is the multi-GB rank.
    //
    // The sparse accumulation below is the SAME arithmetic. The mapping rule is
    // subtle enough to spell out, because the old code got it via brute force
    // (it scattered the full dense block, both triangles):
    //   * a block DIAGONAL entry (row == column) lands once on the mapped
    //     diagonal;
    //   * an off-diagonal entry whose two coordinates MERGE onto one unique
    //     coordinate (mappedRow == mappedColumn) gets BOTH symmetric halves, so
    //     it contributes 2*value to that diagonal;
    //   * otherwise the two halves land on the two symmetric off-diagonal slots,
    //     which are the same slot in upper-triangle storage, so it contributes
    //     value once.
    // Structural zeros are simply never created; entries that accumulate to zero
    // are kept, so the pattern is the union of the contributions and does not
    // depend on numerical cancellation.
    std::vector<std::vector<std::pair<int, double>>> accumulated(
        static_cast<std::size_t>(dimension));
    for (std::size_t blockIndex = 0; blockIndex < blocks.size(); ++blockIndex) {
        const SymmetricCSR &block = *blocks[blockIndex];
        if (maps[blockIndex].size() != static_cast<std::size_t>(block.dimension)) {
            message = "compatibility map does not cover a reduced block";
            return -2;
        }
        const std::vector<int> &map = maps[blockIndex];
        for (int row = 0; row < block.dimension; ++row) {
            for (int entry = block.rowOffsets[static_cast<std::size_t>(row)];
                 entry < block.rowOffsets[static_cast<std::size_t>(row + 1)];
                 ++entry) {
                const int column =
                    block.columnIndices[static_cast<std::size_t>(entry)];
                const double value =
                    block.upperValues[static_cast<std::size_t>(entry)];
                const int mappedRow = map[static_cast<std::size_t>(row)];
                const int mappedColumn = map[static_cast<std::size_t>(column)];
                if (mappedRow < 0 || mappedRow >= dimension ||
                    mappedColumn < 0 || mappedColumn >= dimension) {
                    message = "compatibility map points outside the assembled "
                              "dimension";
                    return -3;
                }
                if (row == column) {
                    accumulated[static_cast<std::size_t>(mappedRow)]
                        .emplace_back(mappedRow, value);
                } else if (mappedRow == mappedColumn) {
                    accumulated[static_cast<std::size_t>(mappedRow)]
                        .emplace_back(mappedRow, 2.0 * value);
                } else {
                    const int lower = std::min(mappedRow, mappedColumn);
                    const int upper = std::max(mappedRow, mappedColumn);
                    accumulated[static_cast<std::size_t>(lower)]
                        .emplace_back(upper, value);
                }
            }
        }
    }

    result = SymmetricCSR{};
    result.dimension = dimension;
    result.rowOffsets.push_back(0);
    for (int row = 0; row < dimension; ++row) {
        std::vector<std::pair<int, double>> &entries =
            accumulated[static_cast<std::size_t>(row)];
        std::sort(entries.begin(), entries.end(),
                  [](const std::pair<int, double> &left,
                     const std::pair<int, double> &right) {
                      return left.first < right.first;
                  });
        int currentColumn = -1;
        for (const std::pair<int, double> &item : entries) {
            if (item.first != currentColumn) {
                currentColumn = item.first;
                result.columnIndices.push_back(currentColumn);
                result.upperValues.push_back(item.second);
            } else {
                result.upperValues.back() += item.second;
            }
        }
        // Release the row as we go: the accumulator never has to coexist with
        // the finished CSR at full size.
        std::vector<std::pair<int, double>>().swap(entries);
        if (result.columnIndices.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            message = "assembled hierarchy matrix exceeds CSR integer range";
            return -4;
        }
        result.rowOffsets.push_back(static_cast<int>(result.columnIndices.size()));
    }
    for (const double value : result.upperValues) {
        if (!std::isfinite(value)) {
            message = "dense hierarchy matrix contains a non-finite value";
            return -5;
        }
    }
    return result.validate(message);
}

int denseGeneralizedSolve(
    const SymmetricCSR &stiffness,
    const SymmetricCSR &mass,
    int numberOfModes,
    std::vector<double> &eigenvalues,
    std::vector<double> &eigenvectors,
    std::string &message)
{
    if (stiffness.dimension < 1 || mass.dimension != stiffness.dimension ||
        numberOfModes < 1 || numberOfModes > stiffness.dimension) {
        message = "invalid final dense generalized eigenproblem dimensions";
        return -1;
    }
    DenseMatrix stiffnessDense = denseFromCSR(stiffness);
    DenseMatrix massDense = denseFromCSR(mass);
    int order = stiffness.dimension;
    int itype = 1;
    char jobz = 'V';
    char range = 'I';
    char uplo = 'U';
    int leading = order;
    double lowerValue = 0.0;
    double upperValue = 0.0;
    int first = 1;
    int last = numberOfModes;
    double absoluteTolerance = 0.0;
    int found = 0;
    eigenvalues.assign(static_cast<std::size_t>(order), 0.0);
    eigenvectors.assign(static_cast<std::size_t>(order) * numberOfModes, 0.0);
    int leadingVectors = order;
    int workspaceSize = std::max(1, 8 * order);
    std::vector<double> workspace(static_cast<std::size_t>(workspaceSize));
    std::vector<int> integerWorkspace(static_cast<std::size_t>(5 * order));
    std::vector<int> failed(static_cast<std::size_t>(order));
    int info = 0;
    LADRUNO_HIERARCHY_DSYGVX(
        &itype, &jobz, &range, &uplo, &order,
        stiffnessDense.values.data(), &leading,
        massDense.values.data(), &leading,
        &lowerValue, &upperValue, &first, &last, &absoluteTolerance,
        &found, eigenvalues.data(), eigenvectors.data(), &leadingVectors,
        workspace.data(), &workspaceSize, integerWorkspace.data(), failed.data(), &info);
    if (info != 0 || found != numberOfModes) {
        message = "final LAPACK dsygvx failed with INFO=" + std::to_string(info) +
            " found=" + std::to_string(found);
        return -2;
    }
    eigenvalues.resize(static_cast<std::size_t>(found));
    return 0;
}

double vectorNorm(const std::vector<double> &values)
{
    double square = 0.0;
    for (double value : values)
        square += value * value;
    return std::sqrt(std::max(0.0, square));
}

double pencilResidual(
    const SymmetricCSR &stiffness,
    const SymmetricCSR &mass,
    const std::vector<double> &vector,
    double eigenvalue)
{
    const std::vector<double> stiffnessVector = stiffness.multiply(vector);
    const std::vector<double> massVector = mass.multiply(vector);
    std::vector<double> residual(vector.size());
    for (std::size_t index = 0; index < vector.size(); ++index)
        residual[index] = stiffnessVector[index] - eigenvalue * massVector[index];
    return vectorNorm(residual) /
        std::max(std::numeric_limits<double>::min(),
                 vectorNorm(stiffnessVector) + std::fabs(eigenvalue) * vectorNorm(massVector));
}

bool closeValue(double left, double right, double relativeTolerance, double absoluteTolerance)
{
    return std::fabs(left - right) <= absoluteTolerance +
        relativeTolerance * std::max(1.0, std::max(std::fabs(left), std::fabs(right)));
}

using SparseEntries = std::map<std::pair<int, int>, double>;

int csrFromSparseEntries(
    int dimension,
    const SparseEntries &entries,
    SymmetricCSR &result,
    std::string &message)
{
    result = SymmetricCSR{};
    if (dimension < 0 ||
        entries.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        message = "sparse entry dimensions exceed the supported integer range";
        return -1;
    }
    result.dimension = dimension;
    result.rowOffsets.assign(static_cast<std::size_t>(dimension) + 1u, 0);
    for (const auto &entry : entries) {
        const int row = entry.first.first;
        const int column = entry.first.second;
        if (row < 0 || row > column || column >= dimension ||
            !std::isfinite(entry.second)) {
            message = "invalid sparse entry while constructing CSR";
            return -2;
        }
        if (entry.second != 0.0)
            ++result.rowOffsets[static_cast<std::size_t>(row + 1)];
    }
    for (int row = 0; row < dimension; ++row)
        result.rowOffsets[static_cast<std::size_t>(row + 1)] +=
            result.rowOffsets[static_cast<std::size_t>(row)];
    result.columnIndices.reserve(entries.size());
    result.upperValues.reserve(entries.size());
    for (const auto &entry : entries) {
        if (entry.second == 0.0)
            continue;
        result.columnIndices.push_back(entry.first.second);
        result.upperValues.push_back(entry.second);
    }
    return result.validate(message);
}

void addCSRToEntries(
    const SymmetricCSR &matrix,
    const std::vector<int> *equations,
    SparseEntries &entries)
{
    for (int row = 0; row < matrix.dimension; ++row) {
        const int globalRow = equations == nullptr
            ? row
            : (*equations)[static_cast<std::size_t>(row)];
        for (int position = matrix.rowOffsets[static_cast<std::size_t>(row)];
             position < matrix.rowOffsets[static_cast<std::size_t>(row + 1)]; ++position) {
            const int column = matrix.columnIndices[static_cast<std::size_t>(position)];
            const int globalColumn = equations == nullptr
                ? column
                : (*equations)[static_cast<std::size_t>(column)];
            entries[{std::min(globalRow, globalColumn),
                     std::max(globalRow, globalColumn)}] +=
                matrix.upperValues[static_cast<std::size_t>(position)];
        }
    }
}

bool sameSparseValues(
    const SparseEntries &assembled,
    const SymmetricCSR &reference,
    double relativeTolerance,
    double absoluteTolerance)
{
    SparseEntries referenceEntries;
    addCSRToEntries(reference, nullptr, referenceEntries);
    for (const auto &entry : assembled) {
        const auto found = referenceEntries.find(entry.first);
        const double expected = found == referenceEntries.end() ? 0.0 : found->second;
        if (!closeValue(entry.second, expected, relativeTolerance, absoluteTolerance))
            return false;
    }
    for (const auto &entry : referenceEntries) {
        const auto found = assembled.find(entry.first);
        const double actual = found == assembled.end() ? 0.0 : found->second;
        if (!closeValue(actual, entry.second, relativeTolerance, absoluteTolerance))
            return false;
    }
    return true;
}

int backSubstitute(
    const DenseMatrix &finalCoordinates,
    int globalDimension,
    const std::vector<FineReduction> &fineReductions,
    const std::vector<Level1Assembly> &level1Assemblies,
    const std::vector<CoarseReduction> &coarseReductions,
    const std::vector<std::vector<int>> &coarseToFinal,
    DenseMatrix &globalCoordinates,
    double &maximumDuplicateJump,
    std::string &message)
{
    globalCoordinates = makeDense(globalDimension, finalCoordinates.columns);
    std::vector<int> contributionCounts(static_cast<std::size_t>(globalDimension), 0);
    for (std::size_t coarseIndex = 0; coarseIndex < coarseReductions.size(); ++coarseIndex) {
        const CoarseReduction &coarse = coarseReductions[coarseIndex];
        const Level1Assembly &assembly =
            level1Assemblies[static_cast<std::size_t>(coarse.assemblyIndex)];
        DenseMatrix coarseReduced = makeDense(
            coarse.cb.stiffness.dimension, finalCoordinates.columns);
        for (int row = 0; row < coarseReduced.rows; ++row) {
            const int finalRow = coarseToFinal[coarseIndex][static_cast<std::size_t>(row)];
            for (int column = 0; column < coarseReduced.columns; ++column)
                coarseReduced(row, column) = finalCoordinates(finalRow, column);
        }
        const DenseMatrix level1Coordinates =
            multiply(coarse.cb.transformation, coarseReduced);

        for (std::size_t childPosition = 0; childPosition < assembly.fineIndices.size();
             ++childPosition) {
            const FineReduction &fine =
                fineReductions[static_cast<std::size_t>(assembly.fineIndices[childPosition])];
            DenseMatrix fineReduced = makeDense(
                fine.cb.stiffness.dimension, finalCoordinates.columns);
            for (int row = 0; row < fineReduced.rows; ++row) {
                const int level1Row =
                    assembly.fineToUnique[childPosition][static_cast<std::size_t>(row)];
                for (int column = 0; column < fineReduced.columns; ++column)
                    fineReduced(row, column) = level1Coordinates(level1Row, column);
            }
            const DenseMatrix localPhysical = multiply(fine.cb.transformation, fineReduced);
            for (int localRow = 0; localRow < localPhysical.rows; ++localRow) {
                const int equation = fine.input->equations[static_cast<std::size_t>(localRow)];
                const int previousCount = contributionCounts[static_cast<std::size_t>(equation)];
                for (int column = 0; column < localPhysical.columns; ++column) {
                    const double value = localPhysical(localRow, column);
                    if (previousCount > 0) {
                        maximumDuplicateJump = std::max(
                            maximumDuplicateJump,
                            std::fabs(value - globalCoordinates(equation, column) /
                                static_cast<double>(previousCount)));
                    }
                    globalCoordinates(equation, column) += value;
                }
                ++contributionCounts[static_cast<std::size_t>(equation)];
            }
        }
    }
    for (int equation = 0; equation < globalDimension; ++equation) {
        const int count = contributionCounts[static_cast<std::size_t>(equation)];
        if (count < 1) {
            message = "hierarchical back-substitution missed a global equation";
            return -1;
        }
        for (int column = 0; column < globalCoordinates.columns; ++column)
            globalCoordinates(equation, column) /= static_cast<double>(count);
    }
    return 0;
}

int checkedTriangularCount(int dimension, int &count)
{
    if (dimension < 0)
        return -1;
    const long long value =
        static_cast<long long>(dimension) * (static_cast<long long>(dimension) + 1) / 2;
    if (value > INT_MAX)
        return -1;
    count = static_cast<int>(value);
    return 0;
}

int makeDisplacements(
    const std::vector<int> &counts,
    std::vector<int> &displacements,
    int &total)
{
    displacements.assign(counts.size(), 0);
    long long running = 0;
    for (std::size_t index = 0; index < counts.size(); ++index) {
        if (counts[index] < 0 || running > INT_MAX) {
            total = 0;
            return -1;
        }
        displacements[index] = static_cast<int>(running);
        running += counts[index];
    }
    if (running > INT_MAX) {
        total = 0;
        return -1;
    }
    total = static_cast<int>(running);
    return 0;
}

int gatherCompatiblePencil(
    const std::vector<CoordinateKey> &localKeys,
    const SymmetricCSR &localStiffness,
    const SymmetricCSR &localMass,
    MPI_Comm communicator,
    std::vector<int> &localToUnique,
    std::vector<std::vector<int>> &mapsOnRoot,
    std::vector<CoordinateKey> &uniqueOnRoot,
    std::vector<int> &countsOnRoot,
    SymmetricCSR &stiffnessOnRoot,
    SymmetricCSR &massOnRoot,
    std::string &message)
{
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);
    const int localDimension = localStiffness.dimension;
    int localFailure = localDimension < 1 || localMass.dimension != localDimension ||
        localKeys.size() != static_cast<std::size_t>(localDimension) ||
        localStiffness.validate(message) < 0 || localMass.validate(message) < 0;
    int globalFailure = 0;
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX, communicator);
    if (globalFailure) {
        message = "a compatibility participant supplied an invalid reduced block";
        return -1;
    }

    std::vector<int> dimensions(rank == 0 ? static_cast<std::size_t>(size) : 0u);
    MPI_Gather(&localDimension, 1, MPI_INT,
               rank == 0 ? dimensions.data() : nullptr, 1, MPI_INT, 0, communicator);
    int localTriangular = 0;
    if (checkedTriangularCount(localDimension, localTriangular) < 0)
        localFailure = 1;
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX, communicator);
    if (globalFailure) {
        message = "a reduced block exceeds MPI triangular count range";
        return -2;
    }

    std::vector<int> encodedKeys(static_cast<std::size_t>(3 * localDimension));
    for (int row = 0; row < localDimension; ++row) {
        const CoordinateKey &key = localKeys[static_cast<std::size_t>(row)];
        encodedKeys[static_cast<std::size_t>(3 * row)] = static_cast<int>(key.kind);
        encodedKeys[static_cast<std::size_t>(3 * row + 1)] = key.owner;
        encodedKeys[static_cast<std::size_t>(3 * row + 2)] = key.index;
    }
    const DenseMatrix denseK = denseFromCSR(localStiffness);
    const DenseMatrix denseM = denseFromCSR(localMass);
    std::vector<double> packedK(static_cast<std::size_t>(localTriangular));
    std::vector<double> packedM(static_cast<std::size_t>(localTriangular));
    int packed = 0;
    for (int column = 0; column < localDimension; ++column) {
        for (int row = 0; row <= column; ++row, ++packed) {
            packedK[static_cast<std::size_t>(packed)] = denseK(row, column);
            packedM[static_cast<std::size_t>(packed)] = denseM(row, column);
        }
    }

    std::vector<int> keyCounts;
    std::vector<int> keyDisplacements;
    std::vector<int> valueCounts;
    std::vector<int> valueDisplacements;
    std::vector<int> mapDisplacements;
    int totalKeys = 0;
    int totalValues = 0;
    if (rank == 0) {
        keyCounts.resize(static_cast<std::size_t>(size));
        valueCounts.resize(static_cast<std::size_t>(size));
        for (int participant = 0; participant < size; ++participant) {
            const long long keyCount = 3LL * dimensions[static_cast<std::size_t>(participant)];
            if (keyCount > INT_MAX || checkedTriangularCount(
                    dimensions[static_cast<std::size_t>(participant)],
                    valueCounts[static_cast<std::size_t>(participant)]) < 0) {
                localFailure = 1;
                break;
            }
            keyCounts[static_cast<std::size_t>(participant)] = static_cast<int>(keyCount);
        }
        if (!localFailure &&
            (makeDisplacements(keyCounts, keyDisplacements, totalKeys) < 0 ||
             makeDisplacements(valueCounts, valueDisplacements, totalValues) < 0))
            localFailure = 1;
        int totalMaps = 0;
        if (!localFailure &&
            makeDisplacements(dimensions, mapDisplacements, totalMaps) < 0)
            localFailure = 1;
    }
    MPI_Bcast(&localFailure, 1, MPI_INT, 0, communicator);
    if (localFailure) {
        message = "compatibility gather counts exceed MPI int range";
        return -3;
    }
    std::vector<int> gatheredKeys(rank == 0 ? static_cast<std::size_t>(totalKeys) : 0u);
    std::vector<double> gatheredK(rank == 0 ? static_cast<std::size_t>(totalValues) : 0u);
    std::vector<double> gatheredM(rank == 0 ? static_cast<std::size_t>(totalValues) : 0u);
    MPI_Gatherv(encodedKeys.data(), 3 * localDimension, MPI_INT,
                rank == 0 ? gatheredKeys.data() : nullptr,
                rank == 0 ? keyCounts.data() : nullptr,
                rank == 0 ? keyDisplacements.data() : nullptr,
                MPI_INT, 0, communicator);
    MPI_Gatherv(packedK.data(), localTriangular, MPI_DOUBLE,
                rank == 0 ? gatheredK.data() : nullptr,
                rank == 0 ? valueCounts.data() : nullptr,
                rank == 0 ? valueDisplacements.data() : nullptr,
                MPI_DOUBLE, 0, communicator);
    MPI_Gatherv(packedM.data(), localTriangular, MPI_DOUBLE,
                rank == 0 ? gatheredM.data() : nullptr,
                rank == 0 ? valueCounts.data() : nullptr,
                rank == 0 ? valueDisplacements.data() : nullptr,
                MPI_DOUBLE, 0, communicator);

    mapsOnRoot.clear();
    uniqueOnRoot.clear();
    countsOnRoot.clear();
    localFailure = 0;
    if (rank == 0) {
        std::vector<std::vector<CoordinateKey>> keyGroups;
        std::vector<SymmetricCSR> stiffnessBlocks;
        std::vector<SymmetricCSR> massBlocks;
        keyGroups.reserve(static_cast<std::size_t>(size));
        stiffnessBlocks.reserve(static_cast<std::size_t>(size));
        massBlocks.reserve(static_cast<std::size_t>(size));
        for (int participant = 0; participant < size; ++participant) {
            const int dimension = dimensions[static_cast<std::size_t>(participant)];
            std::vector<CoordinateKey> keys;
            keys.reserve(static_cast<std::size_t>(dimension));
            const int keyOffset = keyDisplacements[static_cast<std::size_t>(participant)];
            for (int row = 0; row < dimension; ++row) {
                const int kind = gatheredKeys[static_cast<std::size_t>(keyOffset + 3 * row)];
                if (kind < static_cast<int>(KeyKind::Level2Mode) ||
                    kind > static_cast<int>(KeyKind::Level1Mode)) {
                    localFailure = 1;
                    break;
                }
                keys.push_back({
                    static_cast<KeyKind>(kind),
                    gatheredKeys[static_cast<std::size_t>(keyOffset + 3 * row + 1)],
                    gatheredKeys[static_cast<std::size_t>(keyOffset + 3 * row + 2)]});
            }
            DenseMatrix blockK = makeDense(dimension, dimension);
            DenseMatrix blockM = makeDense(dimension, dimension);
            int entry = valueDisplacements[static_cast<std::size_t>(participant)];
            for (int column = 0; column < dimension; ++column) {
                for (int row = 0; row <= column; ++row, ++entry) {
                    blockK(row, column) = blockK(column, row) =
                        gatheredK[static_cast<std::size_t>(entry)];
                    blockM(row, column) = blockM(column, row) =
                        gatheredM[static_cast<std::size_t>(entry)];
                }
            }
            SymmetricCSR blockKCSR;
            SymmetricCSR blockMCSR;
            if (localFailure || csrFromDenseUpper(blockK, blockKCSR, message) < 0 ||
                csrFromDenseUpper(blockM, blockMCSR, message) < 0) {
                localFailure = 1;
                break;
            }
            keyGroups.push_back(std::move(keys));
            stiffnessBlocks.push_back(std::move(blockKCSR));
            massBlocks.push_back(std::move(blockMCSR));
        }
        if (!localFailure && compatibilityMaps(
                keyGroups, mapsOnRoot, uniqueOnRoot, countsOnRoot) < 0)
            localFailure = 1;
        if (!localFailure) {
            std::vector<const SymmetricCSR *> stiffnessPointers;
            std::vector<const SymmetricCSR *> massPointers;
            for (int participant = 0; participant < size; ++participant) {
                stiffnessPointers.push_back(&stiffnessBlocks[static_cast<std::size_t>(participant)]);
                massPointers.push_back(&massBlocks[static_cast<std::size_t>(participant)]);
            }
            if (assembleCompatible(
                    stiffnessPointers, mapsOnRoot,
                    static_cast<int>(uniqueOnRoot.size()), stiffnessOnRoot, message) < 0 ||
                assembleCompatible(
                    massPointers, mapsOnRoot,
                    static_cast<int>(uniqueOnRoot.size()), massOnRoot, message) < 0)
                localFailure = 1;
        }
    }
    MPI_Bcast(&localFailure, 1, MPI_INT, 0, communicator);
    if (localFailure) {
        message = "compatibility assembly failed on its leader";
        return -4;
    }

    std::vector<int> flatMaps;
    if (rank == 0) {
        flatMaps.reserve(static_cast<std::size_t>(totalKeys / 3));
        for (const auto &map : mapsOnRoot)
            flatMaps.insert(flatMaps.end(), map.begin(), map.end());
    }
    localToUnique.resize(static_cast<std::size_t>(localDimension));
    MPI_Scatterv(rank == 0 ? flatMaps.data() : nullptr,
                 rank == 0 ? dimensions.data() : nullptr,
                 rank == 0 ? mapDisplacements.data() : nullptr,
                 MPI_INT, localToUnique.data(), localDimension, MPI_INT, 0, communicator);
    return 0;
}

} // namespace

int solveTwoLevelHierarchy(
    const TwoLevelHierarchyInput &input,
    TwoLevelHierarchyResult &result,
    std::string &message)
{
    message.clear();
    result = TwoLevelHierarchyResult{};
    const int globalDimension = input.globalStiffness.dimension;
    if (globalDimension < 1 || input.globalMass.dimension != globalDimension ||
        input.globalStiffness.validate(message) < 0 || input.globalMass.validate(message) < 0 ||
        input.fineSubdomains.empty() || input.numberOfCoarseSubdomains < 1 ||
        input.numberOfModes < 1 || input.denseMax < input.numberOfModes ||
        !(input.localTolerance > 0.0)) {
        if (message.empty())
            message = "invalid two-level hierarchy dimensions or controls";
        return -1;
    }
    result.diagnostics.originalDimension = globalDimension;

    std::vector<std::set<int>> fineOwners(static_cast<std::size_t>(globalDimension));
    std::vector<std::set<int>> coarseOwners(static_cast<std::size_t>(globalDimension));
    std::set<int> fineIDs;
    for (std::size_t inputIndex = 0; inputIndex < input.fineSubdomains.size(); ++inputIndex) {
        const FineSubdomainInput &fine = input.fineSubdomains[inputIndex];
        if (fine.fine < 0 || fine.coarse < 0 || fine.coarse >= input.numberOfCoarseSubdomains ||
            fine.equations.empty() || fine.stiffness.dimension != static_cast<int>(fine.equations.size()) ||
            fine.mass.dimension != fine.stiffness.dimension || !fineIDs.insert(fine.fine).second ||
            fine.stiffness.validate(message) < 0 || fine.mass.validate(message) < 0) {
            if (!message.empty())
                message = "invalid fine-subdomain matrix: " + message;
            else
                message = "invalid or duplicate fine-subdomain descriptor";
            return -2;
        }
        std::set<int> localEquations;
        for (int equation : fine.equations) {
            if (equation < 0 || equation >= globalDimension ||
                !localEquations.insert(equation).second) {
                message = "fine subdomain contains an invalid or duplicate global equation";
                return -3;
            }
            fineOwners[static_cast<std::size_t>(equation)].insert(fine.fine);
            coarseOwners[static_cast<std::size_t>(equation)].insert(fine.coarse);
        }
    }
    for (int equation = 0; equation < globalDimension; ++equation) {
        if (fineOwners[static_cast<std::size_t>(equation)].empty()) {
            message = "at least one global equation is absent from all fine subdomains";
            return -4;
        }
    }

    // Verify that disjoint local contributions reproduce the supplied global pencil.
    SparseEntries assembledStiffness;
    SparseEntries assembledMass;
    for (const FineSubdomainInput &fine : input.fineSubdomains) {
        addCSRToEntries(fine.stiffness, &fine.equations, assembledStiffness);
        addCSRToEntries(fine.mass, &fine.equations, assembledMass);
    }
    if (!sameSparseValues(assembledStiffness, input.globalStiffness,
                          input.assemblyRtol, input.assemblyAtol) ||
        !sameSparseValues(assembledMass, input.globalMass,
                          input.assemblyRtol, input.assemblyAtol)) {
        message = "fine-subdomain contributions do not reproduce the global pencil";
        return -5;
    }

    std::vector<FineReduction> fineReductions;
    fineReductions.reserve(input.fineSubdomains.size());
    for (const FineSubdomainInput &fine : input.fineSubdomains) {
        std::vector<int> boundary;
        std::vector<int> interior;
        for (std::size_t local = 0; local < fine.equations.size(); ++local) {
            const int equation = fine.equations[local];
            if (fineOwners[static_cast<std::size_t>(equation)].size() > 1u)
                boundary.push_back(static_cast<int>(local));
            else
                interior.push_back(static_cast<int>(local));
        }
        FineReduction reduction;
        reduction.input = &fine;
        if (reduceCraigBampton(
                fine.stiffness, fine.mass, interior, boundary, fine.modesToKeep,
                input, 2, fine.fine, reduction.cb, message) < 0)
            return -6;
        for (int mode = 0; mode < reduction.cb.retainedModes; ++mode)
            reduction.keys.push_back({KeyKind::Level2Mode, fine.fine, mode});
        for (int localBoundary : boundary)
            reduction.keys.push_back(
                {KeyKind::PhysicalEquation, -1,
                 fine.equations[static_cast<std::size_t>(localBoundary)]});
        result.diagnostics.afterLevel2BeforeCompatibility += reduction.cb.stiffness.dimension;
        fineReductions.push_back(std::move(reduction));
    }
    result.diagnostics.appliedT2 = true;

    std::vector<Level1Assembly> level1Assemblies;
    for (int coarse = 0; coarse < input.numberOfCoarseSubdomains; ++coarse) {
        Level1Assembly assembly;
        assembly.coarse = coarse;
        std::vector<std::vector<CoordinateKey>> keyGroups;
        std::vector<const SymmetricCSR *> stiffnessBlocks;
        std::vector<const SymmetricCSR *> massBlocks;
        for (std::size_t fineIndex = 0; fineIndex < fineReductions.size(); ++fineIndex) {
            if (fineReductions[fineIndex].input->coarse != coarse)
                continue;
            assembly.fineIndices.push_back(static_cast<int>(fineIndex));
            keyGroups.push_back(fineReductions[fineIndex].keys);
            stiffnessBlocks.push_back(&fineReductions[fineIndex].cb.stiffness);
            massBlocks.push_back(&fineReductions[fineIndex].cb.mass);
        }
        if (assembly.fineIndices.empty()) {
            message = "a coarse subdomain has no fine children";
            return -7;
        }
        if (compatibilityMaps(
                keyGroups, assembly.fineToUnique, assembly.keys, assembly.counts) < 0)
            return -8;
        if (assembleCompatible(
                stiffnessBlocks, assembly.fineToUnique,
                static_cast<int>(assembly.keys.size()), assembly.stiffness, message) < 0 ||
            assembleCompatible(
                massBlocks, assembly.fineToUnique,
                static_cast<int>(assembly.keys.size()), assembly.mass, message) < 0)
            return -9;
        result.diagnostics.afterLevel2Compatibility += assembly.stiffness.dimension;
        result.diagnostics.level2CompatibilityCounts.insert(
            result.diagnostics.level2CompatibilityCounts.end(),
            assembly.counts.begin(), assembly.counts.end());
        level1Assemblies.push_back(std::move(assembly));
    }
    result.diagnostics.appliedS2 = true;

    std::vector<CoarseReduction> coarseReductions;
    for (std::size_t assemblyIndex = 0; assemblyIndex < level1Assemblies.size(); ++assemblyIndex) {
        const Level1Assembly &assembly = level1Assemblies[assemblyIndex];
        std::vector<int> boundary;
        std::vector<int> interior;
        for (std::size_t coordinate = 0; coordinate < assembly.keys.size(); ++coordinate) {
            const CoordinateKey &key = assembly.keys[coordinate];
            if (key.kind == KeyKind::PhysicalEquation &&
                coarseOwners[static_cast<std::size_t>(key.index)].size() > 1u)
                boundary.push_back(static_cast<int>(coordinate));
            else
                interior.push_back(static_cast<int>(coordinate));
        }
        CoarseReduction reduction;
        reduction.coarse = assembly.coarse;
        reduction.assemblyIndex = static_cast<int>(assemblyIndex);
        if (input.diagnosticAblateLevel1) {
            // BENCHMARK DIAGNOSTIC ONLY (ADR-1000 P3d / P4 section 2b). Skip the
            // level-1 Craig-Bampton reduction and hand S1 the un-reduced group
            // pencil, which yields the "level-2 only" space. This is NOT a
            // solver configuration and NOT a fallback: no command option can set
            // this flag, and LadrunoCMSEigenSolver refuses an input that carries
            // it. The reconstruction path needs cb.transformation, which for an
            // un-reduced group is exactly the identity.
            reduction.cb.stiffness = assembly.stiffness;
            reduction.cb.mass = assembly.mass;
            reduction.cb.interior = interior;
            reduction.cb.boundary = boundary;
            reduction.cb.retainedModes = 0;
            reduction.cb.transformation =
                makeDense(assembly.stiffness.dimension, assembly.stiffness.dimension);
            if (reduction.cb.transformation.values.empty() &&
                assembly.stiffness.dimension > 0) {
                message = "level-1 ablation could not allocate its identity "
                          "transformation";
                return -10;
            }
            for (int coordinate = 0; coordinate < assembly.stiffness.dimension;
                 ++coordinate)
                reduction.cb.transformation(coordinate, coordinate) = 1.0;
            reduction.keys = assembly.keys;
        } else {
            if (reduceCraigBampton(
                    assembly.stiffness, assembly.mass, interior, boundary,
                    input.modesLevel1, input, 1, assembly.coarse,
                    reduction.cb, message) < 0)
                return -10;
            for (int mode = 0; mode < reduction.cb.retainedModes; ++mode)
                reduction.keys.push_back({KeyKind::Level1Mode, assembly.coarse, mode});
            for (int boundaryCoordinate : boundary)
                reduction.keys.push_back(assembly.keys[static_cast<std::size_t>(boundaryCoordinate)]);
        }
        result.diagnostics.afterLevel1BeforeCompatibility += reduction.cb.stiffness.dimension;
        coarseReductions.push_back(std::move(reduction));
    }
    // Labelled, not silently equivalent: an ablated run reports appliedT1 =
    // false, so no downstream consumer can mistake it for the mandatory chain.
    result.diagnostics.appliedT1 = !input.diagnosticAblateLevel1;
    result.diagnostics.ablatedLevel1 = input.diagnosticAblateLevel1;

    std::vector<std::vector<CoordinateKey>> coarseKeyGroups;
    std::vector<const SymmetricCSR *> coarseStiffnessBlocks;
    std::vector<const SymmetricCSR *> coarseMassBlocks;
    for (const CoarseReduction &coarse : coarseReductions) {
        coarseKeyGroups.push_back(coarse.keys);
        coarseStiffnessBlocks.push_back(&coarse.cb.stiffness);
        coarseMassBlocks.push_back(&coarse.cb.mass);
    }
    std::vector<std::vector<int>> coarseToFinal;
    std::vector<CoordinateKey> finalKeys;
    if (compatibilityMaps(
            coarseKeyGroups, coarseToFinal, finalKeys,
            result.diagnostics.level1CompatibilityCounts) < 0)
        return -11;
    if (assembleCompatible(
            coarseStiffnessBlocks, coarseToFinal, static_cast<int>(finalKeys.size()),
            result.finalStiffness, message) < 0 ||
        assembleCompatible(
            coarseMassBlocks, coarseToFinal, static_cast<int>(finalKeys.size()),
            result.finalMass, message) < 0)
        return -12;
    result.diagnostics.appliedS1 = true;
    result.diagnostics.finalRawDimension = result.finalStiffness.dimension;

    StaticCondensation finalCondensation;
    if (condenseCoordinateMass(
            result.finalStiffness, result.finalMass,
            input.massRtol, input.massAtol,
            finalCondensation, message) < 0) {
        message = "final hierarchy mass condensation failed: " + message;
        return -13;
    }
    result.diagnostics.finalMasslessDimension =
        static_cast<int>(finalCondensation.massless.size());
    result.diagnostics.finalDynamicDimension =
        static_cast<int>(finalCondensation.dynamic.size());
    if (result.diagnostics.finalDynamicDimension > input.denseMax) {
        message = "final dynamic dimension exceeds denseMax";
        return -14;
    }
    if (input.numberOfModes > result.diagnostics.finalDynamicDimension) {
        message = "final dynamic hierarchy dimension is smaller than numberOfModes";
        return -15;
    }
    std::vector<double> dynamicVectors;
    if (denseGeneralizedSolve(
            finalCondensation.reducedStiffness,
            finalCondensation.reducedMass,
            input.numberOfModes,
            result.eigenvalues,
            dynamicVectors,
            message) < 0)
        return -16;
    std::vector<double> finalVectors;
    if (reconstructStaticCoordinates(
            finalCondensation, dynamicVectors, input.numberOfModes,
            finalVectors, message) < 0)
        return -17;

    const int finalDimension = result.finalStiffness.dimension;
    DenseMatrix finalEigenvectors;
    finalEigenvectors.rows = finalDimension;
    finalEigenvectors.columns = input.numberOfModes;
    finalEigenvectors.values = std::move(finalVectors);
    DenseMatrix globalEigenvectors;
    if (backSubstitute(
            finalEigenvectors, globalDimension, fineReductions, level1Assemblies,
            coarseReductions, coarseToFinal, globalEigenvectors,
            result.diagnostics.maximumDuplicateJump, message) < 0)
        return -18;
    result.eigenvectors = globalEigenvectors.values;

    if (input.storeTotalTransformation) {
        std::size_t transformationEntries = 0;
        if (!checkedProduct(static_cast<std::size_t>(globalDimension),
                            static_cast<std::size_t>(finalDimension),
                            transformationEntries) ||
            transformationEntries > input.maximumTransformationEntries) {
            message = "requested diagnostic transformation exceeds its entry limit";
            return -19;
        }
        DenseMatrix finalIdentity = makeDense(finalDimension, finalDimension);
        for (int diagonal = 0; diagonal < finalDimension; ++diagonal)
            finalIdentity(diagonal, diagonal) = 1.0;
        DenseMatrix totalTransformation;
        double diagnosticJump = 0.0;
        if (backSubstitute(
                finalIdentity, globalDimension, fineReductions, level1Assemblies,
                coarseReductions, coarseToFinal, totalTransformation,
                diagnosticJump, message) < 0)
            return -20;
        result.totalTransformation = std::move(totalTransformation.values);
        result.diagnostics.maximumDuplicateJump = std::max(
            result.diagnostics.maximumDuplicateJump, diagnosticJump);
    }
    result.normalizedResiduals.assign(static_cast<std::size_t>(input.numberOfModes), 0.0);
    result.diagnostics.maximumResidualEquation.assign(
        static_cast<std::size_t>(input.numberOfModes), -1);
    result.diagnostics.maximumResidualRole.assign(
        static_cast<std::size_t>(input.numberOfModes), -1);
    result.diagnostics.maximumResidualMassless.assign(
        static_cast<std::size_t>(input.numberOfModes), -1);
    result.diagnostics.residualNorm.assign(
        static_cast<std::size_t>(input.numberOfModes), 0.0);
    result.diagnostics.stiffnessActionNorm.assign(
        static_cast<std::size_t>(input.numberOfModes), 0.0);
    result.diagnostics.massActionNorm.assign(
        static_cast<std::size_t>(input.numberOfModes), 0.0);
    for (int mode = 0; mode < input.numberOfModes; ++mode) {
        std::vector<double> vector(static_cast<std::size_t>(globalDimension));
        for (int row = 0; row < globalDimension; ++row)
            vector[static_cast<std::size_t>(row)] = globalEigenvectors(row, mode);
        result.normalizedResiduals[static_cast<std::size_t>(mode)] = pencilResidual(
            input.globalStiffness, input.globalMass, vector,
            result.eigenvalues[static_cast<std::size_t>(mode)]);
    }
    return 0;
}

int solveDistributedHierarchy(
    const DistributedHierarchyInput &input,
    DistributedHierarchyResult &result,
    std::string &message)
{
    result = DistributedHierarchyResult{};
    message.clear();
    int rank = 0;
    int worldSize = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    int localFailure = input.globalDimension < 1 || input.fine != rank ||
        input.coarse < 0 || input.coarse >= input.numberOfCoarseSubdomains ||
        input.localEquations.empty() ||
        input.localStiffness.dimension != static_cast<int>(input.localEquations.size()) ||
        input.localMass.dimension != input.localStiffness.dimension ||
        input.ownedStiffnessContributions == nullptr ||
        input.ownedMassContributions == nullptr ||
        input.modesLevel2 < 1 || input.modesLevel1 < 1 || input.numberOfModes < 1 ||
        input.denseMax < input.numberOfModes || !(input.tolerance > 0.0);
    int globalFailure = 0;
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (globalFailure) {
        message = "invalid distributed hierarchy input on at least one rank";
        return -1;
    }

    // Cross-rank option-consistency gate. numberOfModes, denseMax, modesLevel2
    // and modesLevel1 size the numberOfModes-wide collectives below (the
    // MPI_Bcast of the eigenvalues, the per-mode reductions). If ranks disagree
    // — a heterogeneous launch with divergent option files — those collectives
    // are undefined behavior (count mismatch => silent corruption or hang). The
    // per-rank validation above only checks each rank locally, so verify
    // agreement collectively (symmetric on every rank) before any such
    // collective runs. This is additive and cannot itself deadlock.
    {
        const int consistencyValues[4] = {
            input.numberOfModes, input.denseMax,
            input.modesLevel2, input.modesLevel1};
        int consistencyMin[4];
        int consistencyMax[4];
        MPI_Allreduce(consistencyValues, consistencyMin, 4, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
        MPI_Allreduce(consistencyValues, consistencyMax, 4, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD);
        for (int i = 0; i < 4; ++i) {
            if (consistencyMin[i] != consistencyMax[i]) {
                static const char *const names[4] = {
                    "numberOfModes", "denseMax", "modesLevel2", "modesLevel1"};
                message = std::string("distributed hierarchy option '") +
                    names[i] + "' disagrees across ranks (min=" +
                    std::to_string(consistencyMin[i]) + ", max=" +
                    std::to_string(consistencyMax[i]) +
                    "); every rank must request the same value";
                return -1;
            }
        }
    }

    MPI_Comm groupComm = MPI_COMM_NULL;
    MPI_Comm leaderComm = MPI_COMM_NULL;
    MPI_Comm_split(MPI_COMM_WORLD, input.coarse, rank, &groupComm);
    int groupRank = 0;
    int groupSize = 0;
    MPI_Comm_rank(groupComm, &groupRank);
    MPI_Comm_size(groupComm, &groupSize);
    MPI_Comm_split(MPI_COMM_WORLD, groupRank == 0 ? 0 : MPI_UNDEFINED,
                   rank, &leaderComm);
    const auto cleanup = [&]() {
        if (leaderComm != MPI_COMM_NULL)
            MPI_Comm_free(&leaderComm);
        if (groupComm != MPI_COMM_NULL)
            MPI_Comm_free(&groupComm);
    };

    TwoLevelHierarchyInput controls;
    // The physical hierarchy is followed by original-pencil refinement, so
    // forcing every fixed-interface mode to be ten times tighter than the
    // requested final tolerance can reject an otherwise adequate CMS basis
    // at a Lanczos roundoff plateau.  Keep the local solve at least as strict
    // as the public tolerance (and never looser than 1e-8); the final
    // residual gate remains authoritative.
    controls.localTolerance = std::min(1.0e-8, input.tolerance);
    controls.maximumOperatorApplications = input.maximumOperatorApplications;
    controls.maximumRestarts = 20;
    controls.massRtol = input.massRtol;
    controls.massAtol = input.massAtol;

    std::vector<int> localSeen(static_cast<std::size_t>(input.globalDimension), 0);
    localFailure = 0;
    for (int equation : input.localEquations) {
        if (equation < 0 || equation >= input.globalDimension ||
            localSeen[static_cast<std::size_t>(equation)] != 0) {
            localFailure = 1;
            break;
        }
        localSeen[static_cast<std::size_t>(equation)] = 1;
    }
    std::vector<int> fineOwnerCounts(static_cast<std::size_t>(input.globalDimension), 0);
    MPI_Allreduce(localSeen.data(), fineOwnerCounts.data(), input.globalDimension,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (std::find(fineOwnerCounts.begin(), fineOwnerCounts.end(), 0) !=
        fineOwnerCounts.end())
        localFailure = 1;

    std::vector<int> groupSeen(static_cast<std::size_t>(input.globalDimension), 0);
    MPI_Allreduce(localSeen.data(), groupSeen.data(), input.globalDimension,
                  MPI_INT, MPI_MAX, groupComm);
    std::vector<int> leaderSeen(static_cast<std::size_t>(input.globalDimension), 0);
    if (groupRank == 0)
        leaderSeen = groupSeen;
    std::vector<int> coarseOwnerCounts(static_cast<std::size_t>(input.globalDimension), 0);
    MPI_Allreduce(leaderSeen.data(), coarseOwnerCounts.data(), input.globalDimension,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (globalFailure) {
        cleanup();
        message = "distributed local equation ownership is incomplete or invalid";
        return -2;
    }

    // ADR-1000 P4 section 1: per-phase attribution of the distributed
    // hierarchy. Each phase stamp is taken AFTER the collective that closes the
    // phase, so a phase's time includes the wait for the slowest rank in it --
    // which is what you want when hunting a bottleneck across ranks.
    double phaseMark = MPI_Wtime();
    const auto stampPhase = [&phaseMark]() {
        const double now = MPI_Wtime();
        const double elapsed = now - phaseMark;
        phaseMark = now;
        return elapsed;
    };

    std::vector<int> interior;
    std::vector<int> boundary;
    for (std::size_t local = 0; local < input.localEquations.size(); ++local) {
        const int equation = input.localEquations[local];
        if (fineOwnerCounts[static_cast<std::size_t>(equation)] > 1)
            boundary.push_back(static_cast<int>(local));
        else
            interior.push_back(static_cast<int>(local));
    }
    result.diagnostics.partitionSeconds = stampPhase();

    CBReduction fineReduction;
    localFailure = reduceCraigBampton(
        input.localStiffness, input.localMass, interior, boundary,
        input.modesLevel2, controls, 2, input.fine,
        fineReduction, message) < 0;
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (globalFailure) {
        cleanup();
        if (!localFailure)
            message = "level-2 reduction failed on another rank";
        return -3;
    }
    result.diagnostics.appliedT2 = true;
    result.diagnostics.fineModesSeconds = stampPhase();
    result.diagnostics.originalDimension = input.globalDimension;
    int localDimension = fineReduction.stiffness.dimension;
    MPI_Allreduce(&localDimension,
                  &result.diagnostics.afterLevel2BeforeCompatibility,
                  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    std::vector<CoordinateKey> fineKeys;
    for (int mode = 0; mode < fineReduction.retainedModes; ++mode)
        fineKeys.push_back({KeyKind::Level2Mode, input.fine, mode});
    for (int localBoundary : boundary)
        fineKeys.push_back({
            KeyKind::PhysicalEquation, -1,
            input.localEquations[static_cast<std::size_t>(localBoundary)]});

    std::vector<int> fineToGroup;
    std::vector<std::vector<int>> fineMapsOnLeader;
    std::vector<CoordinateKey> groupKeys;
    std::vector<int> groupCounts;
    SymmetricCSR groupStiffness;
    SymmetricCSR groupMass;
    localFailure = gatherCompatiblePencil(
        fineKeys, fineReduction.stiffness, fineReduction.mass, groupComm,
        fineToGroup, fineMapsOnLeader, groupKeys, groupCounts,
        groupStiffness, groupMass, message) < 0;
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (globalFailure) {
        cleanup();
        if (!localFailure)
            message = "level-2 compatibility failed in another coarse group";
        return -4;
    }
    result.diagnostics.appliedS2 = true;
    result.diagnostics.compatibilitySeconds = stampPhase();
    int groupCompatibleDimension = groupRank == 0 ? groupStiffness.dimension : 0;
    MPI_Allreduce(&groupCompatibleDimension,
                  &result.diagnostics.afterLevel2Compatibility,
                  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    CBReduction coarseReduction;
    std::vector<CoordinateKey> coarseKeys;
    localFailure = 0;
    if (groupRank == 0) {
        std::vector<int> coarseInterior;
        std::vector<int> coarseBoundary;
        for (std::size_t coordinate = 0; coordinate < groupKeys.size(); ++coordinate) {
            const CoordinateKey &key = groupKeys[coordinate];
            if (key.kind == KeyKind::PhysicalEquation &&
                coarseOwnerCounts[static_cast<std::size_t>(key.index)] > 1)
                coarseBoundary.push_back(static_cast<int>(coordinate));
            else
                coarseInterior.push_back(static_cast<int>(coordinate));
        }
        if (reduceCraigBampton(
                groupStiffness, groupMass, coarseInterior, coarseBoundary,
                input.modesLevel1, controls, 1, input.coarse,
                coarseReduction, message) < 0) {
            localFailure = 1;
        } else {
            for (int mode = 0; mode < coarseReduction.retainedModes; ++mode)
                coarseKeys.push_back({KeyKind::Level1Mode, input.coarse, mode});
            for (int coordinate : coarseBoundary)
                coarseKeys.push_back(groupKeys[static_cast<std::size_t>(coordinate)]);
        }
    }
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (globalFailure) {
        cleanup();
        if (!localFailure)
            message = "level-1 reduction failed in another coarse group";
        return -5;
    }
    result.diagnostics.appliedT1 = true;
    result.diagnostics.level1Seconds = stampPhase();
    int coarseDimension = groupRank == 0 ? coarseReduction.stiffness.dimension : 0;
    MPI_Allreduce(&coarseDimension,
                  &result.diagnostics.afterLevel1BeforeCompatibility,
                  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    std::vector<int> coarseToFinal;
    std::vector<std::vector<int>> coarseMapsOnRoot;
    std::vector<CoordinateKey> finalKeys;
    std::vector<int> finalCounts;
    SymmetricCSR finalStiffness;
    SymmetricCSR finalMass;
    localFailure = 0;
    if (groupRank == 0) {
        localFailure = gatherCompatiblePencil(
            coarseKeys, coarseReduction.stiffness, coarseReduction.mass, leaderComm,
            coarseToFinal, coarseMapsOnRoot, finalKeys, finalCounts,
            finalStiffness, finalMass, message) < 0;
    }
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (globalFailure) {
        cleanup();
        if (!localFailure)
            message = "level-1 compatibility failed on the leader communicator";
        return -6;
    }
    result.diagnostics.appliedS1 = true;

    std::vector<double> finalVectors;
    result.eigenvalues.assign(static_cast<std::size_t>(input.numberOfModes), 0.0);
    int finalDimension = 0;
    int finalMassless = 0;
    int finalDynamic = 0;
    localFailure = 0;
    if (rank == 0) {
        finalDimension = finalStiffness.dimension;
        StaticCondensation condensation;
        if (condenseCoordinateMass(
                finalStiffness, finalMass, input.massRtol, input.massAtol,
                condensation, message) < 0) {
            localFailure = 1;
        } else {
            finalMassless = static_cast<int>(condensation.massless.size());
            finalDynamic = static_cast<int>(condensation.dynamic.size());
            if (finalDynamic > input.denseMax || input.numberOfModes > finalDynamic) {
                message = finalDynamic > input.denseMax
                    ? "final dynamic dimension rD=" + std::to_string(finalDynamic) +
                      " exceeds denseMax=" + std::to_string(input.denseMax)
                    : "final dynamic dimension rD=" + std::to_string(finalDynamic) +
                      " is smaller than numberOfModes=" +
                      std::to_string(input.numberOfModes);
                localFailure = 1;
            } else {
                std::vector<double> dynamicVectors;
                if (denseGeneralizedSolve(
                        condensation.reducedStiffness, condensation.reducedMass,
                        input.numberOfModes, result.eigenvalues,
                        dynamicVectors, message) < 0 ||
                    reconstructStaticCoordinates(
                        condensation, dynamicVectors, input.numberOfModes,
                        finalVectors, message) < 0)
                    localFailure = 1;
            }
        }
    }
    MPI_Bcast(&localFailure, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (localFailure) {
        cleanup();
        if (rank != 0)
            message = "final condensed eigensolve failed on rank 0";
        return -7;
    }
    MPI_Bcast(&finalDimension, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&finalMassless, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&finalDynamic, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(result.eigenvalues.data(), input.numberOfModes,
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
    result.diagnostics.finalRawDimension = finalDimension;
    result.diagnostics.finalMasslessDimension = finalMassless;
    result.diagnostics.finalDynamicDimension = finalDynamic;
    int finalVectorCount = 0;
    if (!checkedMPIProduct(
            finalDimension, input.numberOfModes, finalVectorCount)) {
        cleanup();
        message = "final eigenvector broadcast exceeds MPI int range";
        return -8;
    }
    if (groupRank == 0) {
        if (rank != 0)
            finalVectors.resize(static_cast<std::size_t>(finalVectorCount));
        MPI_Bcast(finalVectors.data(), finalVectorCount,
                  MPI_DOUBLE, 0, leaderComm);
    }

    result.diagnostics.globalSolveSeconds = stampPhase();

    std::vector<double> groupCoordinates;
    if (groupRank == 0) {
        DenseMatrix coarseCoordinates = makeDense(
            coarseReduction.stiffness.dimension, input.numberOfModes);
        for (int row = 0; row < coarseCoordinates.rows; ++row) {
            const int finalRow = coarseToFinal[static_cast<std::size_t>(row)];
            for (int mode = 0; mode < input.numberOfModes; ++mode)
                coarseCoordinates(row, mode) = finalVectors[
                    static_cast<std::size_t>(mode) * finalDimension + finalRow];
        }
        groupCoordinates = multiply(
            coarseReduction.transformation, coarseCoordinates).values;
    }

    int receiveCount = 0;
    localFailure = checkedMPIProduct(
        fineReduction.stiffness.dimension, input.numberOfModes,
        receiveCount) ? 0 : 1;
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (globalFailure) {
        cleanup();
        message = "back-substitution receive count exceeds MPI int range";
        return -8;
    }
    std::vector<double> fineCoordinates(static_cast<std::size_t>(receiveCount));
    std::vector<int> sendCounts;
    std::vector<int> sendDisplacements;
    std::vector<double> sendBuffer;
    localFailure = 0;
    if (groupRank == 0) {
        sendCounts.resize(static_cast<std::size_t>(groupSize));
        for (int participant = 0; participant < groupSize; ++participant) {
            const long long count = static_cast<long long>(
                fineMapsOnLeader[static_cast<std::size_t>(participant)].size()) *
                input.numberOfModes;
            if (count > INT_MAX) {
                localFailure = 1;
                break;
            }
            sendCounts[static_cast<std::size_t>(participant)] = static_cast<int>(count);
        }
        int total = 0;
        if (!localFailure &&
            makeDisplacements(sendCounts, sendDisplacements, total) < 0)
            localFailure = 1;
        if (!localFailure) {
            sendBuffer.resize(static_cast<std::size_t>(total));
            const int groupRows = groupStiffness.dimension;
            for (int participant = 0; participant < groupSize; ++participant) {
                const auto &map = fineMapsOnLeader[static_cast<std::size_t>(participant)];
                const int offset = sendDisplacements[static_cast<std::size_t>(participant)];
                for (int mode = 0; mode < input.numberOfModes; ++mode)
                    for (std::size_t row = 0; row < map.size(); ++row)
                        sendBuffer[static_cast<std::size_t>(offset) +
                                   static_cast<std::size_t>(mode) * map.size() + row] =
                            groupCoordinates[
                                static_cast<std::size_t>(mode) * groupRows +
                                map[row]];
            }
        }
    }
    MPI_Bcast(&localFailure, 1, MPI_INT, 0, groupComm);
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (globalFailure) {
        cleanup();
        message = "back-substitution scatter exceeds MPI int range";
        return -8;
    }
    MPI_Scatterv(groupRank == 0 ? sendBuffer.data() : nullptr,
                 groupRank == 0 ? sendCounts.data() : nullptr,
                 groupRank == 0 ? sendDisplacements.data() : nullptr,
                 MPI_DOUBLE, fineCoordinates.data(), receiveCount,
                 MPI_DOUBLE, 0, groupComm);
    DenseMatrix fineReducedCoordinates;
    fineReducedCoordinates.rows = fineReduction.stiffness.dimension;
    fineReducedCoordinates.columns = input.numberOfModes;
    fineReducedCoordinates.values = std::move(fineCoordinates);
    const DenseMatrix localPhysical = multiply(
        fineReduction.transformation, fineReducedCoordinates);

    result.eigenvectors.assign(
        static_cast<std::size_t>(input.globalDimension) * input.numberOfModes, 0.0);
    result.normalizedResiduals.assign(static_cast<std::size_t>(input.numberOfModes), 0.0);
    result.diagnostics.maximumResidualEquation.assign(
        static_cast<std::size_t>(input.numberOfModes), -1);
    result.diagnostics.maximumResidualRole.assign(
        static_cast<std::size_t>(input.numberOfModes), -1);
    result.diagnostics.maximumResidualMassless.assign(
        static_cast<std::size_t>(input.numberOfModes), -1);
    result.diagnostics.residualNorm.assign(
        static_cast<std::size_t>(input.numberOfModes), 0.0);
    result.diagnostics.stiffnessActionNorm.assign(
        static_cast<std::size_t>(input.numberOfModes), 0.0);
    result.diagnostics.massActionNorm.assign(
        static_cast<std::size_t>(input.numberOfModes), 0.0);
    std::vector<double> localValues(static_cast<std::size_t>(input.globalDimension));
    std::vector<double> localSquares(static_cast<std::size_t>(input.globalDimension));
    std::vector<double> sums(static_cast<std::size_t>(input.globalDimension));
    std::vector<double> squareSums(static_cast<std::size_t>(input.globalDimension));
    for (int mode = 0; mode < input.numberOfModes; ++mode) {
        std::fill(localValues.begin(), localValues.end(), 0.0);
        std::fill(localSquares.begin(), localSquares.end(), 0.0);
        for (std::size_t local = 0; local < input.localEquations.size(); ++local) {
            const int equation = input.localEquations[local];
            const double value = localPhysical(static_cast<int>(local), mode);
            localValues[static_cast<std::size_t>(equation)] = value;
            localSquares[static_cast<std::size_t>(equation)] = value * value;
        }
        MPI_Allreduce(localValues.data(), sums.data(), input.globalDimension,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(localSquares.data(), squareSums.data(), input.globalDimension,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (int equation = 0; equation < input.globalDimension; ++equation) {
            const int count = fineOwnerCounts[static_cast<std::size_t>(equation)];
            const double mean = sums[static_cast<std::size_t>(equation)] / count;
            result.eigenvectors[
                static_cast<std::size_t>(mode) * input.globalDimension + equation] = mean;
            const double variance = std::max(
                0.0, squareSums[static_cast<std::size_t>(equation)] / count - mean * mean);
            result.diagnostics.maximumDuplicateJump = std::max(
                result.diagnostics.maximumDuplicateJump, std::sqrt(variance));
        }
    }

    result.diagnostics.backSubstitutionSeconds = stampPhase();

    SymmetricCSR ownedGlobalStiffness;
    SymmetricCSR ownedGlobalMass;
    if (buildSymmetricCSR(
            input.globalDimension, *input.ownedStiffnessContributions,
            ContributionKind::Stiffness, ownedGlobalStiffness, message) < 0 ||
        buildSymmetricCSR(
            input.globalDimension, *input.ownedMassContributions,
            ContributionKind::Mass, ownedGlobalMass, message) < 0) {
        localFailure = 1;
    } else {
        localFailure = 0;
    }
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (globalFailure) {
        cleanup();
        message = "distributed residual pencil assembly failed";
        return -9;
    }
    std::vector<double> globalK(static_cast<std::size_t>(input.globalDimension));
    std::vector<double> globalM(static_cast<std::size_t>(input.globalDimension));
    std::vector<double> localMassDiagonal(static_cast<std::size_t>(input.globalDimension), 0.0);
    std::vector<double> globalMassDiagonal(static_cast<std::size_t>(input.globalDimension), 0.0);
    for (int row = 0; row < ownedGlobalMass.dimension; ++row) {
        for (int position = ownedGlobalMass.rowOffsets[static_cast<std::size_t>(row)];
             position < ownedGlobalMass.rowOffsets[static_cast<std::size_t>(row + 1)];
             ++position) {
            if (ownedGlobalMass.columnIndices[static_cast<std::size_t>(position)] == row)
                localMassDiagonal[static_cast<std::size_t>(row)] +=
                    ownedGlobalMass.upperValues[static_cast<std::size_t>(position)];
        }
    }
    MPI_Allreduce(localMassDiagonal.data(), globalMassDiagonal.data(),
                  input.globalDimension, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double globalMassScale = 1.0;
    for (double value : globalMassDiagonal)
        globalMassScale = std::max(globalMassScale, std::fabs(value));
    const double globalMassThreshold =
        input.massAtol + input.massRtol * globalMassScale;

    // The hierarchical reductions may mix physical D/Z coordinates.  After
    // the complete inverse path, enforce equilibrium once more in the
    // coordinate-aligned nullspace of the original mass matrix.  Otherwise a
    // truncated basis can have accurate Ritz values but very large residuals
    // on massless interface rotations.
    std::vector<int> masslessMap(static_cast<std::size_t>(input.globalDimension), -1);
    std::vector<int> masslessEquations;
    for (int equation = 0; equation < input.globalDimension; ++equation) {
        if (std::fabs(globalMassDiagonal[static_cast<std::size_t>(equation)]) <=
            globalMassThreshold) {
            masslessMap[static_cast<std::size_t>(equation)] =
                static_cast<int>(masslessEquations.size());
            masslessEquations.push_back(equation);
        }
    }
    result.diagnostics.originalMasslessDimension =
        static_cast<int>(masslessEquations.size());
    if (!masslessEquations.empty()) {
        double localMisalignedMass = 0.0;
        for (int row = 0; row < ownedGlobalMass.dimension; ++row) {
            for (int position = ownedGlobalMass.rowOffsets[static_cast<std::size_t>(row)];
                 position < ownedGlobalMass.rowOffsets[static_cast<std::size_t>(row + 1)];
                 ++position) {
                const int column = ownedGlobalMass.columnIndices[
                    static_cast<std::size_t>(position)];
                if (masslessMap[static_cast<std::size_t>(row)] >= 0 ||
                    masslessMap[static_cast<std::size_t>(column)] >= 0)
                    localMisalignedMass = std::max(
                        localMisalignedMass,
                        std::fabs(ownedGlobalMass.upperValues[
                            static_cast<std::size_t>(position)]));
            }
        }
        double globalMisalignedMass = 0.0;
        MPI_Allreduce(&localMisalignedMass, &globalMisalignedMass, 1,
                      MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (globalMisalignedMass > globalMassThreshold) {
            cleanup();
            message = "original mass nullspace is not coordinate aligned";
            return -10;
        }

        SparseEntries localMasslessStiffnessEntries;
        for (int row = 0; row < ownedGlobalStiffness.dimension; ++row) {
            const int masslessRow = masslessMap[static_cast<std::size_t>(row)];
            if (masslessRow < 0)
                continue;
            for (int position = ownedGlobalStiffness.rowOffsets[static_cast<std::size_t>(row)];
                 position < ownedGlobalStiffness.rowOffsets[static_cast<std::size_t>(row + 1)];
                 ++position) {
                const int column = ownedGlobalStiffness.columnIndices[
                    static_cast<std::size_t>(position)];
                const int masslessColumn =
                    masslessMap[static_cast<std::size_t>(column)];
                if (masslessColumn >= 0)
                    localMasslessStiffnessEntries[
                        {std::min(masslessRow, masslessColumn),
                         std::max(masslessRow, masslessColumn)}] +=
                        ownedGlobalStiffness.upperValues[
                            static_cast<std::size_t>(position)];
            }
        }
        SymmetricCSR localMasslessStiffness;
        if (csrFromSparseEntries(
                static_cast<int>(masslessEquations.size()),
                localMasslessStiffnessEntries,
                localMasslessStiffness, message) < 0) {
            cleanup();
            return -11;
        }
        DistributedMumpsSPD masslessFactor;
        if (masslessFactor.factorize(localMasslessStiffness, message) < 0) {
            cleanup();
            message = "original K_ZZ factorization failed: " + message;
            return -12;
        }
        std::vector<double> staticCoordinates(
            masslessEquations.size() * static_cast<std::size_t>(input.numberOfModes),
            0.0);
        std::vector<double> globalDynamicAction(
            static_cast<std::size_t>(input.globalDimension), 0.0);
        for (int mode = 0; mode < input.numberOfModes; ++mode) {
            const auto begin = result.eigenvectors.begin() +
                static_cast<std::ptrdiff_t>(mode) * input.globalDimension;
            std::vector<double> dynamicVector(begin, begin + input.globalDimension);
            for (int equation : masslessEquations)
                dynamicVector[static_cast<std::size_t>(equation)] = 0.0;
            std::vector<double> localDynamicAction =
                ownedGlobalStiffness.multiply(dynamicVector);
            MPI_Allreduce(localDynamicAction.data(), globalDynamicAction.data(),
                          input.globalDimension, MPI_DOUBLE, MPI_SUM,
                          MPI_COMM_WORLD);
            for (std::size_t local = 0; local < masslessEquations.size(); ++local)
                staticCoordinates[
                    static_cast<std::size_t>(mode) * masslessEquations.size() + local] =
                    -globalDynamicAction[
                        static_cast<std::size_t>(masslessEquations[local])];
        }
        if (masslessFactor.solve(
                staticCoordinates, input.numberOfModes, message) < 0) {
            cleanup();
            return -13;
        }
        for (int mode = 0; mode < input.numberOfModes; ++mode) {
            for (std::size_t local = 0; local < masslessEquations.size(); ++local)
                result.eigenvectors[
                    static_cast<std::size_t>(mode) * input.globalDimension +
                    masslessEquations[local]] =
                    staticCoordinates[
                        static_cast<std::size_t>(mode) * masslessEquations.size() + local];
        }

        // Re-establish the Ritz condition in the statically cleaned span.
        const DistributedHierarchyResult cleaned = result;
        if (solveNestedRitzUnion(
                input.globalDimension, input.numberOfModes,
                *input.ownedStiffnessContributions,
                *input.ownedMassContributions,
                cleaned, result, message) < 0) {
            cleanup();
            message = "original D/Z cleanup Ritz solve failed: " + message;
            return -14;
        }

        result.diagnostics.nestedEnrichmentDimension = 0;
    }
    for (int mode = 0; mode < input.numberOfModes; ++mode) {
        const auto begin = result.eigenvectors.begin() +
            static_cast<std::ptrdiff_t>(mode) * input.globalDimension;
        std::vector<double> vector(begin, begin + input.globalDimension);
        const std::vector<double> localK = ownedGlobalStiffness.multiply(vector);
        const std::vector<double> localM = ownedGlobalMass.multiply(vector);
        MPI_Allreduce(localK.data(), globalK.data(), input.globalDimension,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(localM.data(), globalM.data(), input.globalDimension,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        std::vector<double> residual(static_cast<std::size_t>(input.globalDimension));
        for (int equation = 0; equation < input.globalDimension; ++equation)
            residual[static_cast<std::size_t>(equation)] =
                globalK[static_cast<std::size_t>(equation)] -
                result.eigenvalues[static_cast<std::size_t>(mode)] *
                globalM[static_cast<std::size_t>(equation)];
        const auto maximum = std::max_element(
            residual.begin(), residual.end(),
            [](double left, double right) {
                return std::fabs(left) < std::fabs(right);
            });
        const double residualNorm = vectorNorm(residual);
        const double stiffnessNorm = vectorNorm(globalK);
        const double massNorm = vectorNorm(globalM);
        result.diagnostics.maximumResidualEquation[static_cast<std::size_t>(mode)] =
            maximum == residual.end()
                ? -1
                : static_cast<int>(std::distance(residual.begin(), maximum));
        if (maximum != residual.end()) {
            const int equation = static_cast<int>(
                std::distance(residual.begin(), maximum));
            result.diagnostics.maximumResidualRole[static_cast<std::size_t>(mode)] =
                coarseOwnerCounts[static_cast<std::size_t>(equation)] > 1
                    ? 2
                    : (fineOwnerCounts[static_cast<std::size_t>(equation)] > 1 ? 1 : 0);
            result.diagnostics.maximumResidualMassless[static_cast<std::size_t>(mode)] =
                std::fabs(globalMassDiagonal[static_cast<std::size_t>(equation)]) <=
                    globalMassThreshold
                ? 1
                : 0;
        }
        result.diagnostics.residualNorm[static_cast<std::size_t>(mode)] = residualNorm;
        result.diagnostics.stiffnessActionNorm[static_cast<std::size_t>(mode)] = stiffnessNorm;
        result.diagnostics.massActionNorm[static_cast<std::size_t>(mode)] = massNorm;
        result.normalizedResiduals[static_cast<std::size_t>(mode)] =
            residualNorm /
            std::max(std::numeric_limits<double>::min(),
                     stiffnessNorm +
                     std::fabs(result.eigenvalues[static_cast<std::size_t>(mode)]) *
                     massNorm);
    }
    result.diagnostics.publicationSeconds = stampPhase();
    result.diagnostics.peakResidentBytes = peakResidentSetBytes();
    cleanup();
    return 0;
}

int solveNestedRitzUnion(
    int globalDimension,
    int numberOfModes,
    const std::vector<AssemblyRecord> &ownedStiffnessContributions,
    const std::vector<AssemblyRecord> &ownedMassContributions,
    const DistributedHierarchyResult &previous,
    DistributedHierarchyResult &candidate,
    std::string &message)
{
    message.clear();
    const std::size_t vectorEntries =
        static_cast<std::size_t>(globalDimension) * numberOfModes;
    if (globalDimension < 1 || numberOfModes < 1 ||
        previous.eigenvectors.size() != vectorEntries ||
        candidate.eigenvectors.size() != vectorEntries) {
        message = "nested enrichment vector dimensions are inconsistent";
        return -1;
    }

    SymmetricCSR localStiffness;
    SymmetricCSR localMass;
    int localFailure =
        buildSymmetricCSR(globalDimension, ownedStiffnessContributions,
                          ContributionKind::Stiffness, localStiffness, message) < 0 ||
        buildSymmetricCSR(globalDimension, ownedMassContributions,
                          ContributionKind::Mass, localMass, message) < 0;
    int globalFailure = 0;
    MPI_Allreduce(&localFailure, &globalFailure, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    if (globalFailure) {
        if (!localFailure)
            message = "nested enrichment pencil assembly failed on another rank";
        return -2;
    }

    std::vector<std::vector<double>> basis;
    std::vector<std::vector<double>> massBasis;
    basis.reserve(static_cast<std::size_t>(2 * numberOfModes));
    massBasis.reserve(static_cast<std::size_t>(2 * numberOfModes));
    std::vector<double> globalAction(static_cast<std::size_t>(globalDimension));
    const double rankTolerance = 1024.0 * std::numeric_limits<double>::epsilon();
    const DistributedHierarchyResult *sources[2] = {&previous, &candidate};
    for (const DistributedHierarchyResult *source : sources) {
        for (int mode = 0; mode < numberOfModes; ++mode) {
            const auto begin = source->eigenvectors.begin() +
                static_cast<std::ptrdiff_t>(mode) * globalDimension;
            std::vector<double> vector(begin, begin + globalDimension);
            std::vector<double> massVector = localMass.multiply(vector);
            MPI_Allreduce(massVector.data(), globalAction.data(), globalDimension,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            massVector = globalAction;
            double originalMassNorm = 0.0;
            for (int row = 0; row < globalDimension; ++row)
                originalMassNorm += vector[static_cast<std::size_t>(row)] *
                    massVector[static_cast<std::size_t>(row)];

            // Two MGS passes control loss of orthogonality when the two Ritz
            // spaces are nearly coincident.
            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t column = 0; column < basis.size(); ++column) {
                    double projection = 0.0;
                    for (int row = 0; row < globalDimension; ++row)
                        projection += basis[column][static_cast<std::size_t>(row)] *
                            massVector[static_cast<std::size_t>(row)];
                    for (int row = 0; row < globalDimension; ++row) {
                        vector[static_cast<std::size_t>(row)] -=
                            projection * basis[column][static_cast<std::size_t>(row)];
                        massVector[static_cast<std::size_t>(row)] -=
                            projection * massBasis[column][static_cast<std::size_t>(row)];
                    }
                }
            }
            double massNorm = 0.0;
            for (int row = 0; row < globalDimension; ++row)
                massNorm += vector[static_cast<std::size_t>(row)] *
                    massVector[static_cast<std::size_t>(row)];
            if (massNorm <= rankTolerance * std::max(1.0, std::fabs(originalMassNorm)))
                continue;
            const double inverseNorm = 1.0 / std::sqrt(massNorm);
            for (int row = 0; row < globalDimension; ++row) {
                vector[static_cast<std::size_t>(row)] *= inverseNorm;
                massVector[static_cast<std::size_t>(row)] *= inverseNorm;
            }
            basis.push_back(std::move(vector));
            massBasis.push_back(std::move(massVector));
        }
    }
    const int basisDimension = static_cast<int>(basis.size());
    if (basisDimension < numberOfModes) {
        message = "nested enrichment lost rank in the original mass product";
        return -3;
    }

    DenseMatrix projectedStiffness = makeDense(basisDimension, basisDimension);
    std::vector<std::vector<double>> stiffnessBasis(
        static_cast<std::size_t>(basisDimension));
    for (int column = 0; column < basisDimension; ++column) {
        std::vector<double> localAction =
            localStiffness.multiply(basis[static_cast<std::size_t>(column)]);
        MPI_Allreduce(localAction.data(), globalAction.data(), globalDimension,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        stiffnessBasis[static_cast<std::size_t>(column)] = globalAction;
    }
    for (int column = 0; column < basisDimension; ++column) {
        for (int row = 0; row <= column; ++row) {
            double forward = 0.0;
            double reverse = 0.0;
            for (int equation = 0; equation < globalDimension; ++equation) {
                forward += basis[static_cast<std::size_t>(row)]
                    [static_cast<std::size_t>(equation)] *
                    stiffnessBasis[static_cast<std::size_t>(column)]
                    [static_cast<std::size_t>(equation)];
                reverse += basis[static_cast<std::size_t>(column)]
                    [static_cast<std::size_t>(equation)] *
                    stiffnessBasis[static_cast<std::size_t>(row)]
                    [static_cast<std::size_t>(equation)];
            }
            const double value = 0.5 * (forward + reverse);
            projectedStiffness(row, column) = value;
            projectedStiffness(column, row) = value;
        }
    }
    DenseMatrix projectedMass = makeDense(basisDimension, basisDimension);
    for (int diagonal = 0; diagonal < basisDimension; ++diagonal)
        projectedMass(diagonal, diagonal) = 1.0;
    SymmetricCSR stiffnessCSR;
    SymmetricCSR massCSR;
    if (csrFromDenseUpper(projectedStiffness, stiffnessCSR, message) < 0 ||
        csrFromDenseUpper(projectedMass, massCSR, message) < 0)
        return -4;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> coordinates;
    localFailure = 0;
    if (rank == 0 && denseGeneralizedSolve(
            stiffnessCSR, massCSR, numberOfModes,
            candidate.eigenvalues, coordinates, message) < 0)
        localFailure = 1;
    MPI_Bcast(&localFailure, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (localFailure) {
        if (rank != 0)
            message = "nested Rayleigh-Ritz solve failed on rank 0";
        return -5;
    }
    if (rank != 0) {
        candidate.eigenvalues.resize(static_cast<std::size_t>(numberOfModes));
        coordinates.resize(static_cast<std::size_t>(basisDimension) * numberOfModes);
    }
    MPI_Bcast(candidate.eigenvalues.data(), numberOfModes,
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int coordinateCount = 0;
    if (!checkedMPIProduct(basisDimension, numberOfModes, coordinateCount)) {
        message = "nested Rayleigh-Ritz coordinates exceed MPI int range";
        return -6;
    }
    MPI_Bcast(coordinates.data(), coordinateCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    candidate.eigenvectors.assign(vectorEntries, 0.0);
    candidate.normalizedResiduals.assign(
        static_cast<std::size_t>(numberOfModes), 0.0);
    candidate.diagnostics.maximumResidualEquation.assign(
        static_cast<std::size_t>(numberOfModes), -1);
    candidate.diagnostics.maximumResidualRole.assign(
        static_cast<std::size_t>(numberOfModes), -1);
    candidate.diagnostics.maximumResidualMassless.assign(
        static_cast<std::size_t>(numberOfModes), -1);
    candidate.diagnostics.residualNorm.assign(
        static_cast<std::size_t>(numberOfModes), 0.0);
    candidate.diagnostics.stiffnessActionNorm.assign(
        static_cast<std::size_t>(numberOfModes), 0.0);
    candidate.diagnostics.massActionNorm.assign(
        static_cast<std::size_t>(numberOfModes), 0.0);
    for (int mode = 0; mode < numberOfModes; ++mode) {
        std::vector<double> vector(static_cast<std::size_t>(globalDimension), 0.0);
        for (int column = 0; column < basisDimension; ++column) {
            const double coefficient = coordinates[
                static_cast<std::size_t>(mode) * basisDimension + column];
            for (int row = 0; row < globalDimension; ++row)
                vector[static_cast<std::size_t>(row)] += coefficient *
                    basis[static_cast<std::size_t>(column)]
                    [static_cast<std::size_t>(row)];
        }
        std::copy(vector.begin(), vector.end(),
                  candidate.eigenvectors.begin() +
                      static_cast<std::ptrdiff_t>(mode) * globalDimension);
        std::vector<double> localK = localStiffness.multiply(vector);
        std::vector<double> localM = localMass.multiply(vector);
        std::vector<double> globalK(static_cast<std::size_t>(globalDimension));
        std::vector<double> globalM(static_cast<std::size_t>(globalDimension));
        MPI_Allreduce(localK.data(), globalK.data(), globalDimension,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(localM.data(), globalM.data(), globalDimension,
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        std::vector<double> residual(static_cast<std::size_t>(globalDimension));
        for (int row = 0; row < globalDimension; ++row)
            residual[static_cast<std::size_t>(row)] =
                globalK[static_cast<std::size_t>(row)] -
                candidate.eigenvalues[static_cast<std::size_t>(mode)] *
                globalM[static_cast<std::size_t>(row)];
        const auto maximum = std::max_element(
            residual.begin(), residual.end(),
            [](double left, double right) {
                return std::fabs(left) < std::fabs(right);
            });
        const double residualNorm = vectorNorm(residual);
        const double stiffnessNorm = vectorNorm(globalK);
        const double massNorm = vectorNorm(globalM);
        candidate.diagnostics.maximumResidualEquation[static_cast<std::size_t>(mode)] =
            maximum == residual.end()
                ? -1
                : static_cast<int>(std::distance(residual.begin(), maximum));
        candidate.diagnostics.residualNorm[static_cast<std::size_t>(mode)] = residualNorm;
        candidate.diagnostics.stiffnessActionNorm[static_cast<std::size_t>(mode)] = stiffnessNorm;
        candidate.diagnostics.massActionNorm[static_cast<std::size_t>(mode)] = massNorm;
        candidate.normalizedResiduals[static_cast<std::size_t>(mode)] =
            residualNorm /
            std::max(std::numeric_limits<double>::min(),
                     stiffnessNorm +
                     std::fabs(candidate.eigenvalues[static_cast<std::size_t>(mode)]) *
                     massNorm);
    }
    candidate.diagnostics.nestedEnrichmentDimension = basisDimension;
    return 0;
}

} // namespace ladruno_cms
