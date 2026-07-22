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

#ifndef LadrunoCMSEigenSOE_h
#define LadrunoCMSEigenSOE_h

#include "LadrunoCMSOptions.h"

#include <EigenSOE.h>

#include <cstddef>
#include <string>
#include <vector>

class Channel;
class FEM_ObjectBroker;
class LadrunoCMSEigenSolver;

class LadrunoCMSEigenSOE : public EigenSOE
{
public:
    LadrunoCMSEigenSOE(
        LadrunoCMSEigenSolver &solver,
        const ladruno_cms::Options &options);
    ~LadrunoCMSEigenSOE() override = default;

    int getNumEqn() const;
    int setSize(Graph &graph) override;
    int addA(const Matrix &matrix, const ID &id, double factor = 1.0) override;
    int addM(const Matrix &matrix, const ID &id, double factor = 1.0) override;
    void zeroA() override;
    void zeroM() override;

    int sendSelf(int commitTag, Channel &channel) override;
    int recvSelf(int commitTag, Channel &channel, FEM_ObjectBroker &broker) override;

    const ladruno_cms::Options &getOptions() const;
    const std::vector<ladruno_cms::AssemblyRecord> &getStiffnessContributions() const;
    const std::vector<ladruno_cms::AssemblyRecord> &getMassContributions() const;
    const std::vector<int> &getFineLabels() const;
    const std::vector<int> &getCoarseOfFine() const;
    int getWorldRank() const;
    int getWorldSize() const;
    int getNumberOfCoarseSubdomains() const;
    int getFineSubdomainsPerCoarse() const;
    int verifyReplicatedAssembly(std::string &message) const;

private:
    int buildTopology(Graph &graph, std::string &message);
    int addContribution(
        ladruno_cms::ContributionKind kind,
        const Matrix &matrix,
        const ID &id,
        double factor);

    int size_ = 0;
    int worldRank_ = -1;
    int worldSize_ = 0;
    int numberOfCoarseSubdomains_ = 0;
    int fineSubdomainsPerCoarse_ = 0;
    ladruno_cms::Options options_;
    std::vector<int> fineLabels_;
    std::vector<int> coarseOfFine_;
    std::vector<ladruno_cms::AssemblyRecord> stiffness_;
    std::vector<ladruno_cms::AssemblyRecord> mass_;
    std::vector<ladruno_cms::AssemblyRecord> stiffnessSignatures_;
    std::vector<ladruno_cms::AssemblyRecord> massSignatures_;
    std::size_t nextStiffnessOrdinal_ = 0;
    std::size_t nextMassOrdinal_ = 0;
    std::size_t stiffnessVerificationBytes_ = 0;
    std::size_t massVerificationBytes_ = 0;
};

#endif
