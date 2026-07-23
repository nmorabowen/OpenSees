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

#ifndef LadrunoCMSEigenSolver_h
#define LadrunoCMSEigenSolver_h

#include <EigenSolver.h>
#include "LadrunoCMSHierarchy.h"
#include <Vector.h>

#include <vector>

class Channel;
class FEM_ObjectBroker;
class LadrunoCMSEigenSOE;

class LadrunoCMSEigenSolver : public EigenSolver
{
public:
    LadrunoCMSEigenSolver();
    ~LadrunoCMSEigenSolver() override;

    int solve(int numModes, bool generalized, bool findSmallest = true) override;
    int setSize() override;
    int setEigenSOE(LadrunoCMSEigenSOE &soe);
    const Vector &getEigenvector(int mode) override;
    double getEigenvalue(int mode) override;
    const ladruno_cms::HierarchyDiagnostics &getDiagnostics() const;
    const std::vector<double> &getNormalizedResiduals() const;

    int sendSelf(int commitTag, Channel &channel) override;
    int recvSelf(int commitTag, Channel &channel, FEM_ObjectBroker &broker) override;

private:
    LadrunoCMSEigenSOE *soe_ = nullptr;
    std::vector<double> eigenvalues_;
    std::vector<std::vector<double>> eigenvectors_;
    std::vector<double> normalizedResiduals_;
    ladruno_cms::HierarchyDiagnostics diagnostics_;
    Vector *eigenvectorView_ = nullptr;
};

#endif
