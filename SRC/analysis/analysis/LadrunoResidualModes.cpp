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

// Ladruno WP-213: `residualModes` — missing-mass static-correction pseudo-modes
// appended to the domain's mode set (ADR 213 D2/D3).
// See Ladruno_implementation/213_ladruno_modal_basis_completeness_adr.md.
//
//   residualModes -dir d1 <d2 <d3>> <-period rayleigh|$T> <-replace> <-tol $t>
//   residualModes -kinds          (query: 0 = eigen, 1 = residual, 2 = ritz per mode)
//
//   Run after `eigen` (or `eigen -feast`). For each requested global
//   translation d the kernel starts from the rigid shape iota_d (shape start,
//   two blocks), deflated against the current modes:
//       r_d = iota_d - Phi Phi^T M iota_d      (the mass the modes miss)
//       x_d = K^-1 M r_d                       (its static response; K^-1 M Phi
//                                               stays in span(Phi), so this is
//                                               the classical static correction)
//   With r_d in the span the effective modal mass reaches 100 % in direction
//   d; with x_d in the span the static response to M iota_d is reproduced
//   exactly. One static-correction vector alone gives only the second
//   (ADR 213 §3.2: it tops out below 100 % mass). The Rayleigh-Ritz step turns
//   the q <= 6 survivors into M-orthonormal, K-orthogonal pseudo-modes, which
//   are APPENDED to the domain's mode set through the same AnalysisModel calls
//   eigen uses. Every consumer of the mode set (modalProperties,
//   responseSpectrumAnalysis, modalResponseHistory, frequencyResponse,
//   recorders) then uses them with no change; TransformationDOF_Group
//   recovers the constrained nodes.
//
//   Period: every pseudo-mode keeps its Rayleigh-quotient lambda. `-period $T`
//   (ASCE 7-22 §12.9.1.1 single rigid mode: -period 0.05) records T as the
//   mode's Sa-period override in the LadrunoModeKinds registry:
//   responseSpectrumAnalysis then reads Sa(T) for it, while lambda in
//   u = phi Gamma Sa / lambda and the CQC omegas stay the Rayleigh values, so
//   the pseudo-mode contributes the static response times Sa(T) (ADR 213 D3b).

#include <LadrunoRitzKernel.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <EquiSolnAlgo.h>
#include <Domain.h>
#include <Vector.h>
#include <elementAPI.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <cstdlib>
#include <vector>

namespace
{
    const double kTwoPi = 6.283185307179586476925286766559;
}

int
OPS_LadrunoResidualModes(void)
{
    AnalysisModel** modelPtr = OPS_GetAnalysisModel();
    AnalysisModel* theModel = modelPtr ? *modelPtr : 0;
    if (theModel == 0 || theModel->getDomainPtr() == 0) {
        opserr << "residualModes: no AnalysisModel - run eigen first\n";
        return -1;
    }
    Domain* theDomain = theModel->getDomainPtr();

    // A PartitionedDomain / Subdomain (OpenSeesSP) answers false here: its
    // modes and M-products live on the subdomains, which this serial kernel
    // cannot reach (same limitation as modalProperties).
    if (!theDomain->ladrunoThreadedUpdateAllowed()) {
        opserr << "residualModes: a partitioned (OpenSeesSP) domain is not supported\n";
        return -1;
    }

    // ---- parse ----------------------------------------------------------
    std::vector<int> dirs;
    bool replace = false;
    bool fixedPeriod = false;
    double period = 0.0;
    double tol = 1.0e-4;   // a dropped column carries < 1e-8 of the mass (ADR 213 D3)

    // every token is read as text: Python hands integers and floats over as
    // numbers, and OPS_GetStringFromAll is the one reader that accepts both
    // there and leaves classic Tcl's cursor semantics alone
    char tok[64];
    if (OPS_GetNumRemainingInputArgs() == 1) {
        const char* a = OPS_GetStringFromAll(tok, sizeof(tok));
        if (a != 0 && strcmp(a, "-kinds") == 0) {
            std::vector<int> kinds = LadrunoModeKinds::get(theDomain);
            int n = (int)kinds.size();
            if (n > 0 && OPS_SetIntOutput(&n, kinds.data(), false) < 0)
                return -1;
            return 0;
        }
        OPS_ResetCurrentInputArg(-1);
    }

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetStringFromAll(tok, sizeof(tok));
        if (opt == 0) {
            opserr << "residualModes: unreadable argument\n";
            return -1;
        }
        if (strcmp(opt, "-dir") == 0) {
            int numd = 0;
            while (OPS_GetNumRemainingInputArgs() > 0) {
                const char* peek = OPS_GetStringFromAll(tok, sizeof(tok));
                char* end = 0;
                long v = peek ? strtol(peek, &end, 10) : 0;
                if (peek == 0 || end == peek || *end != '\0') {
                    OPS_ResetCurrentInputArg(-1);
                    break;
                }
                const int d = (int)v;
                if (d < 1 || d > 3) {
                    opserr << "residualModes: -dir takes global translations 1, 2 or 3 (got "
                           << d << ")\n";
                    return -1;
                }
                for (size_t k = 0; k < dirs.size(); k++)
                    if (dirs[k] == d) {
                        opserr << "residualModes: direction " << d << " given twice\n";
                        return -1;
                    }
                dirs.push_back(d);
                numd++;
            }
            if (numd == 0) {
                opserr << "residualModes: -dir needs at least one direction\n";
                return -1;
            }
        } else if (strcmp(opt, "-period") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "residualModes: -period needs 'rayleigh' or a period\n";
                return -1;
            }
            const char* p = OPS_GetStringFromAll(tok, sizeof(tok));
            if (p == 0) {
                opserr << "residualModes: -period needs 'rayleigh' or a period\n";
                return -1;
            }
            if (strcmp(p, "rayleigh") == 0) {
                fixedPeriod = false;
            } else {
                char* end = 0;
                double T = strtod(p, &end);
                if (end == p || *end != '\0' || !(T > 0.0) || !std::isfinite(T)) {
                    opserr << "residualModes: -period takes 'rayleigh' or a positive period (got "
                           << p << ")\n";
                    return -1;
                }
                fixedPeriod = true;
                period = T;
            }
        } else if (strcmp(opt, "-replace") == 0) {
            replace = true;
        } else if (strcmp(opt, "-tol") == 0) {
            int one = 1;
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&one, &tol) < 0 ||
                !(tol > 0.0) || !(tol < 1.0)) {
                opserr << "residualModes: -tol needs a value in (0, 1)\n";
                return -1;
            }
        } else {
            opserr << "residualModes: unknown option " << opt << "\n"
                   << "  residualModes -dir d1 <d2 <d3>> <-period rayleigh|$T> <-replace> <-tol $t>\n";
            return -1;
        }
    }
    if (dirs.empty()) {
        opserr << "residualModes: -dir is required\n";
        return -1;
    }

    if (theDomain->getNumEigenvalues() <= 0) {
        opserr << "residualModes: the domain has no modes - run eigen first\n";
        return -1;
    }
    LinearSOE** soePtr = OPS_GetSOE();
    LinearSOE* theSOE = soePtr ? *soePtr : 0;
    if (theSOE == 0) {
        opserr << "residualModes: no LinearSOE - run eigen first\n";
        return -1;
    }

    // ---- current mode set -----------------------------------------------
    LadrunoRitzKernel kernel(*theModel, *theSOE);
    std::vector<int> kinds = LadrunoModeKinds::get(theDomain);
    const Vector& ev = theDomain->getEigenvalues();
    std::vector<double> saPer = LadrunoModeKinds::saPeriods(theDomain);
    std::vector<Vector> phiAll;
    if (kernel.domainModes(phiAll) < 0) {
        opserr << "residualModes: " << kernel.error().c_str() << "\n";
        return -1;
    }

    std::vector<int> masks = LadrunoModeKinds::dirMasks(theDomain);
    std::vector<Vector> phi;
    std::vector<double> lam;
    std::vector<int> keptKinds;
    std::vector<double> keptSaPer;
    std::vector<int> keptMasks;
    int numReplaced = 0;
    int replacedMask = 0;
    for (int k = 0; k < (int)phiAll.size(); k++) {
        if (replace && kinds[k] == LadrunoModeKinds::RESIDUAL) {
            numReplaced++;
            replacedMask |= masks[k];
            continue;
        }
        phi.push_back(phiAll[k]);
        lam.push_back(ev(k));
        keptKinds.push_back(kinds[k]);
        keptSaPer.push_back(saPer[k]);
        keptMasks.push_back(masks[k]);
    }
    if (phi.empty()) {
        opserr << "residualModes: no eigen modes left to complete - run eigen first\n";
        return -1;
    }

    // -replace drops every residual mode, so it rebuilds every direction they
    // covered, not only the ones named in this call
    if (replace) {
        for (int d = 1; d <= 3; d++)
            if ((replacedMask & (1 << (d - 1))) &&
                std::find(dirs.begin(), dirs.end(), d) == dirs.end()) {
                dirs.push_back(d);
                opserr << "residualModes: -replace also rebuilds direction " << d
                       << " (the replaced residual modes covered it)\n";
            }
        std::sort(dirs.begin(), dirs.end());
    }

    // one period rule per residual set: kept residual modes must already use
    // the rule this call asks for, or the deck has to rebuild them (-replace)
    const double wantSa = fixedPeriod ? period : 0.0;
    for (size_t k = 0; k < keptKinds.size(); k++)
        if (keptKinds[k] == LadrunoModeKinds::RESIDUAL && keptSaPer[k] != wantSa) {
            opserr << "residualModes: the set already holds residual modes read at ";
            if (keptSaPer[k] > 0.0)
                opserr << "-period " << keptSaPer[k];
            else
                opserr << "their Rayleigh periods";
            opserr << ", but this call asks for ";
            if (fixedPeriod)
                opserr << "-period " << period;
            else
                opserr << "Rayleigh periods";
            opserr << "; pass -replace to rebuild them all under one rule\n";
            return -1;
        }

    // ---- the kept modes must be eigenpairs of the CURRENT tangent --------
    // A gravity / P-Delta step (or any change of K) after eigen leaves Phi
    // stale, and the pseudo-modes would not be K-orthogonal to it. The Rayleigh
    // quotient of each kept mode on the current K must equal its eigenvalue
    // (exact for the residual modes too: their lambda IS their Rayleigh
    // quotient). Projected on Phi, so Lagrange multiplier rows drop out.
    // Tolerance: an eigenvalue carries round-off of order
    // eps |phi|^T |K| |phi| / phi^T K phi, which Penalty (alpha = 1e14) makes
    // 1e-4 and more on an unchanged model (kip-in, review round 2); the
    // refusal / warning thresholds sit 10x / 1x above that floor, and never
    // below 1e-4 / 1e-7 (WP-214; 100x / 10x let a 20 % stiffness loss after
    // eigen 1 through on the tied-floor model, f1 = 4e-3 vs a measured gap of
    // 2.5e-4). A change smaller than 10x the floor is not refused.
    const int neq = kernel.numEqn();
    {
        std::vector<Vector> Kphi;
        std::vector<double> absQuad;
        kernel.applyK(phi, Kphi, &absQuad);
        Vector Mphi(neq);
        double lamMax = 0.0;
        for (size_t k = 0; k < lam.size(); k++)
            lamMax = std::max(lamMax, std::fabs(lam[k]));
        const double eps = std::numeric_limits<double>::epsilon();
        double worstRefuse = 0.0, worstWarn = 0.0;
        int refuseK = -1, warnK = -1;
        double refuseRel = 0.0, refuseRQ = 0.0, refuseFloor = 0.0, warnRel = 0.0;
        for (size_t k = 0; k < phi.size(); k++) {
            kernel.applyM(phi[k], Mphi);
            const double gm = LadrunoRitzKernel::dot(phi[k], Mphi);
            if (!(gm > 0.0))
                continue;
            const double kq = LadrunoRitzKernel::dot(phi[k], Kphi[k]);
            const double rq = kq / gm;
            const double rel = std::fabs(rq - lam[k]) /
                std::max(std::fabs(lam[k]), 1.0e-12 * lamMax);
            const double floorK = (std::fabs(kq) > 0.0) ? eps * absQuad[k] / std::fabs(kq) : 0.0;
            const double refuseTol = std::max(1.0e-4, 10.0 * floorK);
            const double warnTol = std::max(1.0e-7, floorK);
            if (rel / refuseTol > worstRefuse) {
                worstRefuse = rel / refuseTol;
                refuseK = (int)k;
                refuseRel = rel;
                refuseRQ = rq;
                refuseFloor = floorK;
            }
            if (rel / warnTol > worstWarn) {
                worstWarn = rel / warnTol;
                warnK = (int)k;
                warnRel = rel;
            }
        }
        if (worstRefuse > 1.0) {
            opserr << "residualModes: the model changed after eigen: mode " << refuseK + 1
                   << " has Rayleigh quotient " << refuseRQ << " on the current tangent but"
                   << " eigenvalue " << lam[refuseK] << " (relative " << refuseRel
                   << ", round-off floor " << refuseFloor
                   << "); re-run eigen on the current state first\n";
            return -1;
        }
        if (worstWarn > 1.0)
            opserr << "residualModes: WARNING mode " << warnK + 1
                   << " deviates from the current tangent by " << warnRel
                   << " (relative Rayleigh quotient, above its round-off floor);"
                   << " the eigen solve is loose or the state moved\n";
    }

    // ---- starting vectors M iota_d ---------------------------------------
    std::vector<Vector> F, iota;
    for (size_t k = 0; k < dirs.size(); k++) {
        Vector io(neq), f(neq);
        kernel.rigidInfluence(dirs[k], io);
        kernel.applyM(io, f);
        if (!(f.Norm() > 0.0)) {
            opserr << "residualModes: no free mass in direction " << dirs[k] << "\n";
            return -1;
        }
        iota.push_back(io);
        F.push_back(f);
    }

    // shape starts, two blocks: span{r_d, K^-1 M r_d} with r_d = iota_d
    // deflated against the modes (ADR 213 D3)
    LadrunoRitzKernel::Options opt;
    opt.solveStart = false;
    opt.maxBlocks = 2;
    opt.dropTol = tol;
    LadrunoRitzKernel::Result res;
    const int rc = kernel.run(iota, &phi, opt, res);

    // the kernel overwrote the analysis' LinearSOE with K and its factor:
    // tell the algorithm (WP-188 factor reuse keys on this) before anything
    // can reuse it
    EquiSolnAlgo** algoPtr = OPS_GetAlgorithm();
    if (algoPtr != 0 && *algoPtr != 0)
        (*algoPtr)->domainChanged();

    if (rc < 0) {
        opserr << "residualModes: " << kernel.error().c_str() << "\n";
        return -1;
    }

    const int q = (int)res.psi.size();
    std::vector<double> lamPseudo(q);
    for (int j = 0; j < q; j++)
        lamPseudo[j] = res.lambda[j];   // the mode's own lambda, always

    // ---- report: free mass captured per direction ------------------------
    {
        std::vector<double> gm(phi.size());
        Vector Mphi(neq);
        for (size_t k = 0; k < phi.size(); k++) {
            kernel.applyM(phi[k], Mphi);
            gm[k] = LadrunoRitzKernel::dot(phi[k], Mphi);
        }
        for (size_t d = 0; d < dirs.size(); d++) {
            const double total = LadrunoRitzKernel::dot(iota[d], F[d]);
            double capE = 0.0, capR = 0.0;
            for (size_t k = 0; k < phi.size(); k++)
                if (gm[k] > 0.0) {
                    const double L = LadrunoRitzKernel::dot(phi[k], F[d]);
                    capE += L * L / gm[k];
                }
            for (int j = 0; j < q; j++) {
                const double L = LadrunoRitzKernel::dot(res.psi[j], F[d]);
                capR += L * L;   // psi is M-normalized
            }
            opserr << "residualModes: dir " << dirs[d] << ": " << (int)phi.size()
                   << " modes carry " << 100.0 * capE / total
                   << "% of the free mass, pseudo-modes add " << 100.0 * capR / total
                   << "% -> " << 100.0 * (capE + capR) / total << "%\n";
        }
    }
    if (res.numDropped > 0)
        opserr << "residualModes: " << res.numDropped
               << " residual vector(s) already spanned by the modes (dropped)\n";
    // a kept residual shape this small is close to eigenvector noise (V1:
    // noise sat at ~1e-11 of the mass, the drop rule at 1e-8)
    for (size_t d = 0; d < dirs.size() && d < res.firstBlockRelNorm.size(); d++) {
        const double share = res.firstBlockRelNorm[d] * res.firstBlockRelNorm[d];
        if (share > 0.0 && share < 1.0e-6)
            opserr << "residualModes: WARNING direction " << dirs[d] << ": the modes miss only "
                   << share << " of its mass; the residual shape kept for it may be eigenvector"
                   << " noise (tighten the eigen solve, or raise -tol)\n";
    }

    double lamMaxKept = 0.0;
    for (size_t k = 0; k < lam.size(); k++)
        if (keptKinds[k] == LadrunoModeKinds::EIGEN)
            lamMaxKept = std::max(lamMaxKept, lam[k]);
    for (int j = 0; j < q; j++) {
        const double T = kTwoPi / std::sqrt(lamPseudo[j]);
        opserr << "residualModes: pseudo-mode " << (int)phi.size() + j + 1 << ": T = " << T
               << " s (Rayleigh quotient)";
        if (fixedPeriod)
            opserr << ", RSA reads Sa at T = " << period << " s";
        opserr << "\n";
        if (fixedPeriod && T > period)
            opserr << "residualModes: WARNING pseudo-mode " << (int)phi.size() + j + 1
                   << " is not rigid at " << period << " s (its own period is " << T
                   << " s): Sa(" << period << ") stands in for its real spectral ordinate\n";
        if (res.lambda[j] < lamMaxKept)
            opserr << "residualModes: WARNING pseudo-mode " << (int)phi.size() + j + 1
                   << " has a longer period than the shortest retained mode ("
                   << kTwoPi / std::sqrt(lamMaxKept)
                   << " s): the missing mass is not rigid at this truncation - add modes\n";
    }
    if (fixedPeriod && q > 1)
        opserr << "residualModes: the " << q << " pseudo-modes are ONE rigid mode at " << period
               << " s: responseSpectrumAnalysis -combine sums them before any rule; per-mode"
               << " runs (-mode k, -increment) do not - add their responses algebraically"
               << " before combining with the other modes\n";

    // ---- install the new set ----------------------------------------------
    if (q == 0 && numReplaced == 0) {
        // nothing changes: keep the domain (and its modal properties) as they are
        return 0;
    }
    std::vector<Vector> modes = phi;
    std::vector<double> lambdas = lam;
    std::vector<int> newKinds = keptKinds;
    std::vector<double> newSaPer = keptSaPer;
    std::vector<int> newMasks = keptMasks;
    int callMask = 0;
    for (size_t d = 0; d < dirs.size(); d++)
        callMask |= 1 << (dirs[d] - 1);
    for (int j = 0; j < q; j++) {
        modes.push_back(res.psi[j]);
        lambdas.push_back(lamPseudo[j]);
        newKinds.push_back(LadrunoModeKinds::RESIDUAL);
        newSaPer.push_back(fixedPeriod ? period : 0.0);
        newMasks.push_back(callMask);   // Rayleigh-Ritz mixes the directions
    }
    if (kernel.install(modes, lambdas) < 0) {
        opserr << "residualModes: " << kernel.error().c_str() << "\n";
        return -1;
    }
    LadrunoModeKinds::set(theDomain, newKinds, newSaPer, newMasks);

    // modal properties computed for the old set no longer match it
    theDomain->unsetModalProperties();
    opserr << "residualModes: " << (int)modes.size() << " modes in the set ("
           << q << " residual); re-run modalProperties before RSA\n";

    int nout = q;
    if (q > 0 && OPS_SetDoubleOutput(&nout, lamPseudo.data(), false) < 0)
        return -1;
    return 0;
}
