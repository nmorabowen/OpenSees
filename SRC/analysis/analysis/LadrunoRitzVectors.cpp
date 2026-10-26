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

// Ladruno WP-215: `ritzVectors` — load-dependent Ritz vectors installed as the
// domain's mode set in place of `eigen` (ADR 213 D7).
// See Ladruno_implementation/213_ladruno_modal_basis_completeness_adr.md.
//
//   ritzVectors $nvec <-dir d1 <d2 <d3>>> <-pattern $tag ...> <-cycles $n1 $n2 ...> <-tol $t>
//
//   Starting vectors, one per entry, in the order given on the command line:
//   M iota_d for each -dir d, the reference nodal loads of each -pattern
//   (equation numbering, T^T p under Transformation). The LadrunoRitzKernel
//   runs with load starts (block 0 = K^-1 F) and no basis to deflate against,
//   cycling through the starts block by block (CSI order), until $nvec columns
//   survive; -cycles caps the blocks of each start (CSI `ncyc`; 0 = no cap).
//   The Rayleigh-Ritz pairs (psi, lambda) are installed exactly as eigen
//   installs modes (AnalysisModel::setNumEigenvectors / setEigenvector /
//   setEigenvalues; TransformationDOF_Group recovers the constrained nodes)
//   and recorded as RITZ in LadrunoModeKinds, so modalProperties tags them Z.
//
//   K is the current tangent (K - Kg after a PDelta preload), as eigen sees it.
//   The lowest Ritz values are upper bounds on the lowest eigenvalues; the
//   higher ones are not natural frequencies (WP-215 quirk).
//
//   Report per start f (CARM §19.4.3-4):
//     static  = sum_k (psi_k^T f)^2 / lambda_k  /  f^T K^-1 f   (100 % by construction)
//     dynamic = sum_k (psi_k^T f_m)^2 / f_m^T M^-1 f_m           (f_m: f on the
//               DOFs with mass; for -dir it is the modal mass ratio; for a
//               -pattern it needs a diagonal M and reads n/a otherwise)

#include <LadrunoRitzKernel.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <EquiSolnAlgo.h>
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <Domain.h>
#include <Vector.h>
#include <elementAPI.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
    const double kTwoPi = 6.283185307179586476925286766559;

    struct Start
    {
        bool isDir;   // true: -dir d, false: -pattern tag
        int  value;
    };

    // a whole integer token, or false (the cursor is left on the token)
    bool readIntToken(int& v)
    {
        char tok[64];
        const char* s = OPS_GetStringFromAll(tok, sizeof(tok));
        char* end = 0;
        long l = s ? strtol(s, &end, 10) : 0;
        if (s == 0 || end == s || *end != '\0') {
            OPS_ResetCurrentInputArg(-1);
            return false;
        }
        v = (int)l;
        return true;
    }
}

int
OPS_LadrunoRitzVectors(void)
{
    AnalysisModel** modelPtr = OPS_GetAnalysisModel();
    AnalysisModel* theModel = modelPtr ? *modelPtr : 0;
    StaticAnalysis** saPtr = OPS_GetStaticAnalysis();
    DirectIntegrationAnalysis** taPtr = OPS_GetTransientAnalysis();
    StaticAnalysis* theStatic = saPtr ? *saPtr : 0;
    DirectIntegrationAnalysis* theTransient = taPtr ? *taPtr : 0;
    if (theModel == 0 || theModel->getDomainPtr() == 0 ||
        (theStatic == 0 && theTransient == 0)) {
        opserr << "ritzVectors: no analysis - define one first (e.g. analysis Static,"
               << " with a system that can factor K)\n";
        return -1;
    }
    Domain* theDomain = theModel->getDomainPtr();

    // same limitation as residualModes / modalProperties (OpenSeesSP)
    if (!theDomain->ladrunoThreadedUpdateAllowed()) {
        opserr << "ritzVectors: a partitioned (OpenSeesSP) domain is not supported\n";
        return -1;
    }

    // ---- parse ----------------------------------------------------------
    // every token is read as text (residualModes' rule: Python hands numbers
    // over as numbers, OPS_GetStringFromAll accepts both)
    const char* usage =
        "  ritzVectors $nvec <-dir d1 <d2 <d3>>> <-pattern $tag ...> <-cycles $n1 ...> <-tol $t>\n";
    int nvec = 0;
    if (OPS_GetNumRemainingInputArgs() < 1 || !readIntToken(nvec) || nvec < 1) {
        opserr << "ritzVectors: the number of vectors (>= 1) comes first\n" << usage;
        return -1;
    }
    std::vector<Start> starts;
    std::vector<int> cycles;
    bool haveCycles = false;
    double tol = 1.0e-10;   // the kernel's relative M-norm drop (ADR 213 D1)
    char tok[64];
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetStringFromAll(tok, sizeof(tok));
        if (opt == 0) {
            opserr << "ritzVectors: unreadable argument\n";
            return -1;
        }
        if (strcmp(opt, "-dir") == 0 || strcmp(opt, "-pattern") == 0) {
            const bool isDir = opt[1] == 'd';
            int num = 0, v = 0;
            while (OPS_GetNumRemainingInputArgs() > 0 && readIntToken(v)) {
                if (isDir && (v < 1 || v > 3)) {
                    opserr << "ritzVectors: -dir takes global translations 1, 2 or 3 (got "
                           << v << ")\n";
                    return -1;
                }
                for (size_t k = 0; k < starts.size(); k++)
                    if (starts[k].isDir == isDir && starts[k].value == v) {
                        opserr << "ritzVectors: " << (isDir ? "direction " : "pattern ")
                               << v << " given twice\n";
                        return -1;
                    }
                Start s = {isDir, v};
                starts.push_back(s);
                num++;
            }
            if (num == 0) {
                opserr << "ritzVectors: " << (isDir ? "-dir needs a direction" : "-pattern needs a tag")
                       << "\n";
                return -1;
            }
        } else if (strcmp(opt, "-cycles") == 0) {
            int v = 0;
            haveCycles = true;
            while (OPS_GetNumRemainingInputArgs() > 0 && readIntToken(v)) {
                if (v < 0) {
                    opserr << "ritzVectors: -cycles takes counts >= 0 (0 = no cap)\n";
                    return -1;
                }
                cycles.push_back(v);
            }
        } else if (strcmp(opt, "-tol") == 0) {
            int one = 1;
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&one, &tol) < 0 ||
                !(tol > 0.0) || !(tol < 1.0)) {
                opserr << "ritzVectors: -tol needs a value in (0, 1)\n";
                return -1;
            }
        } else {
            opserr << "ritzVectors: unknown option " << opt << "\n" << usage;
            return -1;
        }
    }
    if (starts.empty()) {
        opserr << "ritzVectors: give at least one starting load (-dir and/or -pattern)\n" << usage;
        return -1;
    }
    if (haveCycles && cycles.size() != starts.size()) {
        opserr << "ritzVectors: -cycles takes one count per starting load, in order ("
               << (int)starts.size() << " starts, " << (int)cycles.size() << " counts)\n";
        return -1;
    }

    // ---- number the model when the domain changed (eigen's rule) ----------
    // ritzVectors stands in for eigen, so it may be the first command after
    // the analysis is defined. Only a moved stamp renumbers: domainChanged()
    // is not free of side effects mid-analysis (CentralDifference resets
    // U(t - dt) to U(t), losing its history; review of #1027)
    int rc = 0;
    const int stamp = theDomain->hasDomainChanged();
    if (theStatic != 0 && theStatic->ladrunoDomainStamp() != stamp)
        rc = theStatic->domainChanged();
    else if (theStatic == 0 && theTransient->ladrunoDomainStamp() != stamp)
        rc = theTransient->domainChanged();
    if (rc < 0) {
        opserr << "ritzVectors: the analysis could not number the model\n";
        return -1;
    }
    LinearSOE** soePtr = OPS_GetSOE();
    LinearSOE* theSOE = soePtr ? *soePtr : 0;
    if (theSOE == 0) {
        opserr << "ritzVectors: the analysis has no LinearSOE\n";
        return -1;
    }

    // ---- starting loads ---------------------------------------------------
    LadrunoRitzKernel kernel(*theModel, *theSOE);
    const int neq = kernel.numEqn();
    std::vector<Vector> F;
    std::vector<std::string> names;
    for (size_t k = 0; k < starts.size(); k++) {
        Vector f(neq);
        char name[64];
        if (starts[k].isDir) {
            Vector io(neq);
            kernel.rigidInfluence(starts[k].value, io);
            kernel.applyM(io, f);
            if (!(f.Norm() > 0.0)) {
                opserr << "ritzVectors: no free mass in direction " << starts[k].value << "\n";
                return -1;
            }
            snprintf(name, sizeof(name), "dir %d", starts[k].value);
        } else {
            int numEle = 0;
            if (kernel.patternLoad(starts[k].value, f, numEle) < 0) {
                opserr << "ritzVectors: " << kernel.error().c_str() << "\n";
                return -1;
            }
            if (numEle > 0)
                opserr << "ritzVectors: WARNING pattern " << starts[k].value << " holds " << numEle
                       << " elemental load(s); only its nodal loads start the vectors\n";
            if (!(f.Norm() > 0.0)) {
                opserr << "ritzVectors: pattern " << starts[k].value
                       << " has no nodal load on a free DOF\n";
                return -1;
            }
            snprintf(name, sizeof(name), "pattern %d", starts[k].value);
        }
        F.push_back(f);
        names.push_back(name);
    }

    // ---- the kernel: load starts, no basis, until nvec columns survive -----
    LadrunoRitzKernel::Options opt;
    opt.solveStart = true;
    opt.maxBlocks = nvec;      // each block keeps >= 1 column or the run ends
    opt.maxVectors = nvec;
    opt.dropTol = tol;
    opt.cyclesPerStart = cycles;
    LadrunoRitzKernel::Result res;
    rc = kernel.run(F, 0, opt, res);

    // the kernel overwrote the analysis' LinearSOE with K and its factor
    // (WP-188 factor reuse keys on this)
    EquiSolnAlgo** algoPtr = OPS_GetAlgorithm();
    if (algoPtr != 0 && *algoPtr != 0)
        (*algoPtr)->domainChanged();

    if (rc < 0) {
        opserr << "ritzVectors: " << kernel.error().c_str() << "\n";
        return -1;
    }
    const int r = (int)res.psi.size();
    if (r == 0) {
        opserr << "ritzVectors: no vector survived (every start is massless?)\n";
        return -1;
    }
    if (r < nvec)
        opserr << "ritzVectors: WARNING only " << r << " of " << nvec << " vectors survived"
               << " (the starts span no more with mass, or -cycles caps them)\n";

    // ---- report: load participation per start ----------------------------
    {
        Vector dM;
        bool offDiag = false;
        kernel.massDiagonal(dM, offDiag);
        for (size_t s = 0; s < F.size(); s++) {
            const Vector& f = F[s];
            Vector x(neq);
            double stat = -1.0;
            if (kernel.solve(f, x) == 0) {
                const double fx = LadrunoRitzKernel::dot(f, x);
                double num = 0.0;
                for (int k = 0; k < r; k++) {
                    const double L = LadrunoRitzKernel::dot(res.psi[k], f);
                    num += L * L / res.lambda[k];
                }
                if (fx > 0.0)
                    stat = num / fx;
            }
            double dyn = -1.0;
            if (starts[s].isDir) {
                Vector io(neq);
                kernel.rigidInfluence(starts[s].value, io);
                const double total = LadrunoRitzKernel::dot(io, f);
                double num = 0.0;
                for (int k = 0; k < r; k++) {
                    const double L = LadrunoRitzKernel::dot(res.psi[k], f);   // psi^T M iota
                    num += L * L;
                }
                if (total > 0.0)
                    dyn = num / total;
            } else if (!offDiag) {
                Vector fm = f;
                double den = 0.0;
                for (int i = 0; i < neq; i++) {
                    if (dM(i) > 0.0)
                        den += fm(i) * fm(i) / dM(i);
                    else
                        fm(i) = 0.0;   // loads on massless DOFs are not counted
                }
                double num = 0.0;
                for (int k = 0; k < r; k++) {
                    const double L = LadrunoRitzKernel::dot(res.psi[k], fm);
                    num += L * L;
                }
                if (den > 0.0)
                    dyn = num / den;
            }
            opserr << "ritzVectors: " << names[s].c_str() << ": " << res.keptPerStart[s]
                   << " vector(s); load participation static ";
            if (stat >= 0.0)
                opserr << 100.0 * stat << "%";
            else
                opserr << "n/a";
            opserr << ", dynamic ";
            if (dyn >= 0.0)
                opserr << 100.0 * dyn << "%";
            else
                opserr << "n/a (consistent mass)";
            opserr << "\n";
        }
    }
    for (int k = 0; k < r; k++)
        opserr << "ritzVectors: vector " << k + 1 << ": T = " << kTwoPi / std::sqrt(res.lambda[k])
               << " s\n";
    opserr << "ritzVectors: only the lowest Ritz periods approximate natural periods (upper"
           << " bounds); do not report the higher ones as natural periods\n";

    // ---- install as the domain's mode set ---------------------------------
    if (kernel.install(res.psi, res.lambda) < 0) {
        opserr << "ritzVectors: " << kernel.error().c_str() << "\n";
        return -1;
    }
    LadrunoModeKinds::set(theDomain, std::vector<int>(r, LadrunoModeKinds::RITZ));
    theDomain->unsetModalProperties();

    int nout = r;
    if (OPS_SetDoubleOutput(&nout, res.lambda.data(), false) < 0)
        return -1;
    return 0;
}
