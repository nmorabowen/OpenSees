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

#include "LadrunoCMSEigenSolver.h"

#include "LadrunoCMSEigenSOE.h"
#include "LadrunoCMSSubspaceRefiner.h"

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <new>
#include <string>

LadrunoCMSEigenSolver::LadrunoCMSEigenSolver()
    : EigenSolver(EigenSOLVER_TAGS_LadrunoCMS)
{
}

LadrunoCMSEigenSolver::~LadrunoCMSEigenSolver()
{
    delete eigenvectorView_;
}

int LadrunoCMSEigenSolver::solve(int numModes, bool generalized, bool findSmallest)
{
    using Clock = std::chrono::steady_clock;
    const auto solveStarted = Clock::now();
    double assemblySeconds = 0.0;
    double hierarchySeconds = 0.0;
    double refinementSeconds = 0.0;
    eigenvalues_.clear();
    eigenvectors_.clear();
    normalizedResiduals_.clear();
    diagnostics_ = ladruno_cms::HierarchyDiagnostics{};
    if (!generalized) {
        opserr << "LadrunoCMSEigenSolver::solve - only generalized problems are supported\n";
        return -1;
    }
    if (!findSmallest) {
        opserr << "LadrunoCMSEigenSolver::solve - findLargest is unsupported\n";
        return -2;
    }
    if (soe_ == nullptr) {
        opserr << "LadrunoCMSEigenSolver::solve - no SOE attached\n";
        return -3;
    }
    int worldSize = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    std::string message;
    const ladruno_cms::Options &options = soe_->getOptions();
    if (options.validate(worldSize, numModes, message) < 0) {
        opserr << "LadrunoCMSEigenSolver::solve - " << message.c_str() << endln;
        return -4;
    }
    if (soe_->verifyReplicatedAssembly(message) < 0) {
        opserr << "LadrunoCMSEigenSolver::solve - " << message.c_str() << endln;
        return -5;
    }
    std::vector<int> localEquations;
    ladruno_cms::SymmetricCSR localStiffness;
    ladruno_cms::SymmetricCSR localMass;
    const auto assemblyStarted = Clock::now();
    if (ladruno_cms::buildOwnedLocalPencil(
            soe_->getStiffnessContributions(), soe_->getMassContributions(),
            localEquations, localStiffness, localMass, message) < 0) {
        opserr << "LadrunoCMSEigenSolver::solve - " << message.c_str() << endln;
        return -6;
    }
    assemblySeconds = std::chrono::duration<double>(
        Clock::now() - assemblyStarted).count();

    const int rank = soe_->getWorldRank();
    const int startingVectors = options.resolvedIterationVectors(numModes, message);
    if (startingVectors < 0) {
        opserr << "LadrunoCMSEigenSolver::solve - " << message.c_str() << endln;
        return -7;
    }
    ladruno_cms::DistributedHierarchyResult accepted;
    bool hierarchySolved = false;
    for (int enrichment = 0; enrichment <= options.maxEnrich; ++enrichment) {
        const long long increment = static_cast<long long>(enrichment) * numModes;
        const long long requestedLevel2 = options.modesLevel2 + increment;
        const long long requestedLevel1 = options.modesLevel1 + increment;
        if (requestedLevel2 > std::numeric_limits<int>::max() ||
            requestedLevel1 > std::numeric_limits<int>::max()) {
            opserr << "LadrunoCMSEigenSolver::solve - enriched mode count "
                   << "exceeds the supported integer range\n";
            return -7;
        }
        ladruno_cms::DistributedHierarchyInput input;
        input.globalDimension = soe_->getNumEqn();
        input.fine = rank;
        input.coarse = soe_->getCoarseOfFine()[static_cast<std::size_t>(rank)];
        input.numberOfCoarseSubdomains = soe_->getNumberOfCoarseSubdomains();
        input.localEquations = localEquations;
        input.localStiffness = localStiffness;
        input.localMass = localMass;
        input.ownedStiffnessContributions = &soe_->getStiffnessContributions();
        input.ownedMassContributions = &soe_->getMassContributions();
        input.modesLevel2 = static_cast<int>(requestedLevel2);
        input.modesLevel1 = static_cast<int>(requestedLevel1);
        input.numberOfModes = startingVectors;
        input.denseMax = options.denseMax;
        input.tolerance = options.tolerance;
        input.maximumOperatorApplications = options.maxIterations;
        input.maximumRestarts = options.maxRestarts;
        input.massRtol = options.massRtol;
        input.massAtol = options.massAtol;
        ladruno_cms::DistributedHierarchyResult candidate;
        const auto hierarchyStarted = Clock::now();
        if (ladruno_cms::solveDistributedHierarchy(input, candidate, message) < 0) {
            hierarchySeconds += std::chrono::duration<double>(
                Clock::now() - hierarchyStarted).count();
            // The "is smaller than numberOfModes" reason is produced on rank 0
            // only; every other rank receives the substitute "failed on rank 0"
            // text, so a rank-local string match desynchronizes the retry
            // decision and the next collective hangs. Broadcast rank 0's
            // verdict so all ranks retry or exit together.
            int insufficientStartingSpace =
                message.find("is smaller than numberOfModes") != std::string::npos
                    ? 1 : 0;
            MPI_Bcast(&insufficientStartingSpace, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (insufficientStartingSpace == 1 && enrichment < options.maxEnrich) {
                if (options.verbose && rank == 0)
                    opserr << "LadrunoCMS: CMS space cannot yet supply q="
                           << startingVectors << "; retrying with k2="
                           << static_cast<int>(requestedLevel2 + numModes)
                           << " k1="
                           << static_cast<int>(requestedLevel1 + numModes)
                           << endln;
                continue;
            }
            opserr << "LadrunoCMSEigenSolver::solve - " << message.c_str() << endln;
            return -8;
        }
        hierarchySeconds += std::chrono::duration<double>(
            Clock::now() - hierarchyStarted).count();
        // ADR-1000 P3d / P4 section 2b: the level-1 ablation is a benchmark-only
        // diagnostic and must never be reachable as a solver configuration or a
        // fallback. No command option can request it and the distributed path
        // never sets it, so this can only fire on a future regression -- which
        // is exactly what it is here to catch. Refuse rather than silently
        // publish modes from a reduction that skipped the mandatory chain.
        if (candidate.diagnostics.ablatedLevel1 ||
            !candidate.diagnostics.appliedT1 || !candidate.diagnostics.appliedS1) {
            opserr << "LadrunoCMSEigenSolver::solve - refusing a result that did "
                      "not apply the mandatory T2 -> S2 -> T1 -> S1 chain "
                      "(level-1 ablation is a benchmark diagnostic, never a "
                      "solver configuration)" << endln;
            return -8;
        }
        const double maximumResidual = *std::max_element(
            candidate.normalizedResiduals.begin(),
            candidate.normalizedResiduals.end());
        if (options.verbose && rank == 0) {
            opserr << "LadrunoCMS initial-space attempt=" << enrichment
                   << " k2=" << static_cast<int>(requestedLevel2)
                   << " k1=" << static_cast<int>(requestedLevel1)
                   << " q=" << startingVectors
                   << " rRaw=" << candidate.diagnostics.finalRawDimension
                   << " rD=" << candidate.diagnostics.finalDynamicDimension
                   << " originalZ=" << candidate.diagnostics.originalMasslessDimension
                   << " maxResidual=" << maximumResidual << endln;
            for (int mode = 0; mode < startingVectors; ++mode) {
                opserr << "  mode=" << mode + 1
                       << " lambda=" << candidate.eigenvalues[static_cast<std::size_t>(mode)]
                       << " residual=" << candidate.normalizedResiduals[static_cast<std::size_t>(mode)]
                       << " maxEq=" << candidate.diagnostics.maximumResidualEquation[static_cast<std::size_t>(mode)]
                       << " role=" << candidate.diagnostics.maximumResidualRole[static_cast<std::size_t>(mode)]
                       << " massless=" << candidate.diagnostics.maximumResidualMassless[static_cast<std::size_t>(mode)]
                       << " ||r||=" << candidate.diagnostics.residualNorm[static_cast<std::size_t>(mode)]
                       << " ||Kx||=" << candidate.diagnostics.stiffnessActionNorm[static_cast<std::size_t>(mode)]
                       << " ||Mx||=" << candidate.diagnostics.massActionNorm[static_cast<std::size_t>(mode)]
                       << endln;
            }
        }
        accepted = std::move(candidate);
        hierarchySolved = true;
        break;
    }
    if (!hierarchySolved) {
        opserr << "LadrunoCMSEigenSolver::solve - CMS hierarchy did not provide "
               << startingVectors << " independent starting vectors\n";
        return -9;
    }

    if (options.refinement == ladruno_cms::RefinementMode::Subspace) {
        const auto refinementStarted = Clock::now();
        ladruno_cms::SymmetricCSR globalStiffness;
        ladruno_cms::SymmetricCSR globalMass;
        if (ladruno_cms::buildSymmetricCSR(
                soe_->getNumEqn(), soe_->getStiffnessContributions(),
                ladruno_cms::ContributionKind::Stiffness,
                globalStiffness, message) < 0 ||
            ladruno_cms::buildSymmetricCSR(
                soe_->getNumEqn(), soe_->getMassContributions(),
                ladruno_cms::ContributionKind::Mass,
                globalMass, message) < 0) {
            opserr << "LadrunoCMSEigenSolver::solve - cannot form the distributed "
                   << "original pencil: " << message.c_str() << endln;
            return -10;
        }
        ladruno_cms::SubspaceRefinementControls controls;
        controls.globalDimension = soe_->getNumEqn();
        controls.numberOfModes = numModes;
        controls.numberOfIterationVectors = startingVectors;
        controls.maximumIterations = options.maxRefineIterations;
        controls.tolerance = options.tolerance;
        controls.massRankTolerance = std::max(
            options.massRtol,
            1024.0 * std::numeric_limits<double>::epsilon());
        ladruno_cms::SubspaceRefinementResult refined;
        const int refineCode = ladruno_cms::refineDistributedSubspace(
            globalStiffness, globalMass, accepted.eigenvectors,
            controls, refined, message);
        if (options.verbose && rank == 0) {
            for (const auto &entry : refined.history)
                opserr << "LadrunoCMS refinement iteration=" << entry.iteration
                       << " maxResidual=" << entry.maximumResidual
                       << " MOrthogonalityLoss=" << entry.massOrthogonalityLoss
                       << endln;
            opserr << "LadrunoCMS refinement timing: factor="
                   << refined.factorizationSeconds
                   << " solve=" << refined.solveSeconds
                   << " inverseRefinement=" << refined.inverseRefinementSeconds
                   << " massAction=" << refined.massActionSeconds
                   << " RayleighRitz=" << refined.rayleighRitzSeconds << endln;
        }
        if (refineCode < 0) {
            opserr << "LadrunoCMSEigenSolver::solve - " << message.c_str() << endln;
            return -11;
        }
        refinementSeconds = std::chrono::duration<double>(
            Clock::now() - refinementStarted).count();
        accepted.eigenvalues = std::move(refined.eigenvalues);
        accepted.eigenvectors = std::move(refined.eigenvectors);
        accepted.normalizedResiduals = std::move(refined.normalizedResiduals);
        accepted.diagnostics.maximumResidualEquation =
            std::move(refined.maximumResidualEquation);
        accepted.diagnostics.residualNorm = std::move(refined.residualNorms);
        accepted.diagnostics.stiffnessActionNorm =
            std::move(refined.stiffnessActionNorms);
        accepted.diagnostics.massActionNorm = std::move(refined.massActionNorms);
        accepted.diagnostics.maximumResidualRole.assign(
            static_cast<std::size_t>(numModes), -1);
        accepted.diagnostics.maximumResidualMassless.assign(
            static_cast<std::size_t>(numModes), -1);
    } else {
        accepted.eigenvalues.resize(static_cast<std::size_t>(numModes));
        accepted.normalizedResiduals.resize(static_cast<std::size_t>(numModes));
        accepted.eigenvectors.resize(
            static_cast<std::size_t>(soe_->getNumEqn()) * numModes);
        accepted.diagnostics.maximumResidualEquation.resize(
            static_cast<std::size_t>(numModes));
        accepted.diagnostics.maximumResidualRole.resize(
            static_cast<std::size_t>(numModes));
        accepted.diagnostics.maximumResidualMassless.resize(
            static_cast<std::size_t>(numModes));
        accepted.diagnostics.residualNorm.resize(static_cast<std::size_t>(numModes));
        accepted.diagnostics.stiffnessActionNorm.resize(
            static_cast<std::size_t>(numModes));
        accepted.diagnostics.massActionNorm.resize(static_cast<std::size_t>(numModes));
        // Uncertified output must warn unconditionally, not only under
        // -verbose: a silent -refine none run is indistinguishable from a
        // certified one at the script level.
        if (rank == 0)
            opserr << "LadrunoCMS: refinement=none publishes an approximate "
                   << "CMS research result without a convergence guarantee\n";
    }
    if (accepted.eigenvalues.size() != static_cast<std::size_t>(numModes) ||
        accepted.eigenvectors.size() !=
            static_cast<std::size_t>(soe_->getNumEqn()) * numModes) {
        opserr << "LadrunoCMSEigenSolver::solve - incomplete accepted result\n";
        return -12;
    }
    eigenvalues_ = std::move(accepted.eigenvalues);
    normalizedResiduals_ = std::move(accepted.normalizedResiduals);
    diagnostics_ = std::move(accepted.diagnostics);
    eigenvectors_.resize(static_cast<std::size_t>(numModes));
    for (int mode = 0; mode < numModes; ++mode) {
        const auto begin = accepted.eigenvectors.begin() +
            static_cast<std::ptrdiff_t>(mode) * soe_->getNumEqn();
        eigenvectors_[static_cast<std::size_t>(mode)].assign(
            begin, begin + soe_->getNumEqn());
    }
    if (options.verbose && rank == 0) {
        opserr << "LadrunoCMS: n=" << diagnostics_.originalDimension
               << " r2=" << diagnostics_.afterLevel2Compatibility
               << " rRaw=" << diagnostics_.finalRawDimension
               << " rD=" << diagnostics_.finalDynamicDimension
               << " maxResidual="
               << *std::max_element(normalizedResiduals_.begin(),
                                    normalizedResiduals_.end()) << endln;
        opserr << "LadrunoCMS timing: assembly=" << assemblySeconds
               << " hierarchy=" << hierarchySeconds
               << " refinement=" << refinementSeconds
               << " total=" << std::chrono::duration<double>(
                      Clock::now() - solveStarted).count() << endln;
        // ADR-1000 P4 section 1: attribute `hierarchy` instead of guessing at
        // it. Phases are rank-local and measured across the collective that
        // closes each one, so they include cross-rank wait -- which is what you
        // want when hunting a bottleneck. Reported on rank 0; use -verbose on a
        // multi-rank run and compare ranks to see imbalance.
        const ladruno_cms::HierarchyDiagnostics &phases = diagnostics_;
        opserr << "LadrunoCMS hierarchy phases [s]: partition="
               << phases.partitionSeconds
               << " T2fineModes=" << phases.fineModesSeconds
               << " S2compatibility=" << phases.compatibilitySeconds
               << " T1level1=" << phases.level1Seconds
               << " S1globalSolve=" << phases.globalSolveSeconds
               << " backSubstitution=" << phases.backSubstitutionSeconds
               << " publication=" << phases.publicationSeconds << endln;
        // Dimension ratios: n/r2 is the level-2 reduction, r2/r* the marginal
        // contribution of level 1 (a ratio near 1.0 means T1 bought nothing).
        const double afterLevel2 =
            static_cast<double>(phases.afterLevel2Compatibility);
        const double finalRaw = static_cast<double>(phases.finalRawDimension);
        opserr << "LadrunoCMS dimensions: n=" << phases.originalDimension
               << " r2=" << phases.afterLevel2Compatibility
               << " rStar=" << phases.finalRawDimension;
        if (afterLevel2 > 0.0)
            opserr << " n/r2="
                   << static_cast<double>(phases.originalDimension) / afterLevel2;
        if (finalRaw > 0.0)
            opserr << " r2/rStar=" << afterLevel2 / finalRaw;
        opserr << endln;
        if (phases.peakResidentBytes > 0u)
            opserr << "LadrunoCMS peak resident set (rank 0): "
                   << static_cast<double>(phases.peakResidentBytes) /
                      (1024.0 * 1024.0) << " MiB" << endln;
    }
    return 0;
}

int LadrunoCMSEigenSolver::setSize()
{
    if (soe_ == nullptr)
        return -1;
    const int size = soe_->getNumEqn();
    if (eigenvectorView_ == nullptr || eigenvectorView_->Size() != size) {
        delete eigenvectorView_;
        eigenvectorView_ = new (std::nothrow) Vector(size);
    }
    eigenvalues_.clear();
    eigenvectors_.clear();
    normalizedResiduals_.clear();
    diagnostics_ = ladruno_cms::HierarchyDiagnostics{};
    return eigenvectorView_ != nullptr && eigenvectorView_->Size() == size ? 0 : -2;
}

int LadrunoCMSEigenSolver::setEigenSOE(LadrunoCMSEigenSOE &soe)
{
    soe_ = &soe;
    return 0;
}

const Vector &LadrunoCMSEigenSolver::getEigenvector(int mode)
{
    if (eigenvectorView_ == nullptr) {
        static Vector empty(0);
        opserr << "LadrunoCMSEigenSolver::getEigenvector - solver is not sized\n";
        return empty;
    }
    eigenvectorView_->Zero();
    if (mode < 1 || static_cast<std::size_t>(mode) > eigenvectors_.size()) {
        opserr << "LadrunoCMSEigenSolver::getEigenvector - mode out of range\n";
        return *eigenvectorView_;
    }
    const std::vector<double> &source = eigenvectors_[static_cast<std::size_t>(mode - 1)];
    if (source.size() != static_cast<std::size_t>(eigenvectorView_->Size())) {
        opserr << "LadrunoCMSEigenSolver::getEigenvector - stored vector size mismatch\n";
        return *eigenvectorView_;
    }
    for (int i = 0; i < eigenvectorView_->Size(); ++i)
        (*eigenvectorView_)(i) = source[static_cast<std::size_t>(i)];
    return *eigenvectorView_;
}

double LadrunoCMSEigenSolver::getEigenvalue(int mode)
{
    if (mode < 1 || static_cast<std::size_t>(mode) > eigenvalues_.size()) {
        opserr << "LadrunoCMSEigenSolver::getEigenvalue - mode out of range\n";
        return 0.0;
    }
    return eigenvalues_[static_cast<std::size_t>(mode - 1)];
}

const ladruno_cms::HierarchyDiagnostics &
LadrunoCMSEigenSolver::getDiagnostics() const
{
    return diagnostics_;
}

const std::vector<double> &
LadrunoCMSEigenSolver::getNormalizedResiduals() const
{
    return normalizedResiduals_;
}

int LadrunoCMSEigenSolver::sendSelf(int, Channel &)
{
    opserr << "LadrunoCMSEigenSolver::sendSelf - unsupported in ADR 1000 v1\n";
    return -1;
}

int LadrunoCMSEigenSolver::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    opserr << "LadrunoCMSEigenSolver::recvSelf - unsupported in ADR 1000 v1\n";
    return -1;
}
