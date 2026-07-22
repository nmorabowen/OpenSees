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

#ifndef LadrunoCMSMumps_h
#define LadrunoCMSMumps_h

#include "LadrunoCMSOptions.h"

#include <memory>
#include <string>
#include <vector>

namespace ladruno_cms {

class MumpsSPD
{
public:
    MumpsSPD();
    ~MumpsSPD();

    MumpsSPD(const MumpsSPD &) = delete;
    MumpsSPD &operator=(const MumpsSPD &) = delete;

    int factorize(const SymmetricCSR &matrix, std::string &message);
    int solve(std::vector<double> &rightHandSides, int numberOfColumns, std::string &message);

    int dimension() const;
    bool isFactorized() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Collective SPD factorization with distributed assembled entries and a
// centralized multi-RHS owned by rank 0.  Every method is collective on
// MPI_COMM_WORLD.
class DistributedMumpsSPD
{
public:
    DistributedMumpsSPD();
    ~DistributedMumpsSPD();

    DistributedMumpsSPD(const DistributedMumpsSPD &) = delete;
    DistributedMumpsSPD &operator=(const DistributedMumpsSPD &) = delete;

    int factorize(const SymmetricCSR &localMatrix, std::string &message);
    int solve(std::vector<double> &rightHandSides, int numberOfColumns,
              std::string &message);
    int dimension() const;
    bool isFactorized() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

struct StaticCondensation {
    StaticCondensation();
    ~StaticCondensation();
    StaticCondensation(StaticCondensation &&) noexcept;
    StaticCondensation &operator=(StaticCondensation &&) noexcept;
    StaticCondensation(const StaticCondensation &) = delete;
    StaticCondensation &operator=(const StaticCondensation &) = delete;

    int originalDimension = 0;
    std::vector<int> dynamic;
    std::vector<int> massless;
    SymmetricCSR reducedStiffness;
    SymmetricCSR reducedMass;

    std::size_t persistentCrossNonzeros() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;

    friend int condenseCoordinateMass(
        const SymmetricCSR &,
        const SymmetricCSR &,
        double,
        double,
        StaticCondensation &,
        std::string &);
    friend int reconstructStaticCoordinates(
        const StaticCondensation &,
        const std::vector<double> &,
        int,
        std::vector<double> &,
        std::string &);
};

int condenseCoordinateMass(
    const SymmetricCSR &stiffness,
    const SymmetricCSR &mass,
    double relativeTolerance,
    double absoluteTolerance,
    StaticCondensation &result,
    std::string &message);

int reconstructStaticCoordinates(
    const StaticCondensation &condensation,
    const std::vector<double> &reducedVectors,
    int numberOfColumns,
    std::vector<double> &fullVectors,
    std::string &message);

} // namespace ladruno_cms

#endif
