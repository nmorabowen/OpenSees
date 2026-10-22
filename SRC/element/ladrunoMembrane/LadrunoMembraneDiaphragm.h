/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

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

#ifndef LadrunoMembraneDiaphragm_h
#define LadrunoMembraneDiaphragm_h

// ============================================================================
// LadrunoMembraneDiaphragm.h — WP-204: in-plane (diaphragm) stiffness of a
// LadrunoMembrane panel, condensed EXACTLY onto every retained node.
//
// The panel is meshed (PSLG CDT from LadrunoMembraneKernel.h), every triangle
// is promoted to a 6-node LST (quadratic, so in-plane bending of a long panel
// converges as h^2 instead of the CST's h), and all non-retained DOFs are
// eliminated by Guyan (static) condensation:
//
//     K_c = K_rr - K_ri K_ii^-1 K_ir,   R = [I; -K_ii^-1 K_ir],   M_c = R' M R
//
// Static condensation is exact for loads applied at the retained nodes, so the
// panel answers like the full fine mesh at its outer vertices, hole vertices,
// interior supports and any extra diaphragm node — the stock condensation kept
// only the outer polygon and slaved every other boundary vertex to straight
// edges (holes, columns and walls inside the panel got zero in-plane
// stiffness). K_ii is solved by a sparse LDL' (Eigen), so a panel with a few
// thousand interior DOFs condenses in milliseconds.
//
// Not included by LadrunoMembraneKernel.h (that header is compiled stand-alone
// by tests/test_ladruno_membrane_kernel_cpp.py, without Eigen on the path).
// ============================================================================

#include "LadrunoMembraneKernel.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace LadrunoMembraneDiaphragm {

using LadrunoMembraneKernel::Point2D;

struct Condensed {
    int nr = 0;                       // retained nodes (2 DOFs each, local e1/e2)
    std::vector<double> K;            // 2nr x 2nr row-major, condensed stiffness
    std::vector<double> M;            // 2nr x 2nr row-major, Guyan mass (rho*t)
    int nTri = 0;
    // strain at each triangle centroid: eps_t = strainOp[t] * u_r,
    // (3 x 2nr row-major per triangle, eps = [exx, eyy, gxy] in e1/e2)
    std::vector<double> strainOp;
    std::vector<double> meanStrainOp; // 3 x 2nr, area-weighted mean over the panel
    std::vector<double> triArea;      // nTri
    std::vector<int> isolated;        // retained slots that touch no triangle
};

namespace detail {

// LST shape-function derivatives with respect to the area coordinates.
// Node order: corners 0,1,2, midsides 3 (0-1), 4 (1-2), 5 (2-0).
inline void lst_dN_dL(const double L[3], double dN[6][3])
{
    for (int i = 0; i < 6; ++i)
        for (int k = 0; k < 3; ++k)
            dN[i][k] = 0.0;
    for (int i = 0; i < 3; ++i)
        dN[i][i] = 4.0 * L[i] - 1.0;
    dN[3][0] = 4.0 * L[1]; dN[3][1] = 4.0 * L[0];
    dN[4][1] = 4.0 * L[2]; dN[4][2] = 4.0 * L[1];
    dN[5][2] = 4.0 * L[0]; dN[5][0] = 4.0 * L[2];
}

inline void lst_N(const double L[3], double N[6])
{
    for (int i = 0; i < 3; ++i)
        N[i] = L[i] * (2.0 * L[i] - 1.0);
    N[3] = 4.0 * L[0] * L[1];
    N[4] = 4.0 * L[1] * L[2];
    N[5] = 4.0 * L[2] * L[0];
}

// B (3 x 12) at area coordinates L; b, c the usual CST gradients * 2A.
inline void lst_B(const double L[3], const double b[3], const double c[3], double A2,
                  double B[3][12])
{
    double dN[6][3];
    lst_dN_dL(L, dN);
    for (int i = 0; i < 6; ++i) {
        double dx = 0.0, dy = 0.0;
        for (int k = 0; k < 3; ++k) {
            dx += dN[i][k] * b[k];
            dy += dN[i][k] * c[k];
        }
        dx /= A2;
        dy /= A2;
        B[0][2 * i] = dx;  B[0][2 * i + 1] = 0.0;
        B[1][2 * i] = 0.0; B[1][2 * i + 1] = dy;
        B[2][2 * i] = dy;  B[2][2 * i + 1] = dx;
    }
}

} // namespace detail

// vertices/triangles: the panel mesh (ccw triangles). retained[r] is the mesh
// vertex of retained slot r (-1: none — the slot gets no stiffness). D is the
// 3x3 plane-stress modulus (row-major, modifiers already applied), t the
// thickness, rhoT the mass per unit area (0: no mass).
inline Condensed condense_lst(const std::vector<Point2D>& vertices,
                              const std::vector<std::array<int, 3>>& triangles,
                              const std::vector<int>& retained,
                              const double D[9], double t, double rhoT)
{
    typedef Eigen::Triplet<double> Trip;
    const int nv = static_cast<int>(vertices.size());
    const int nt = static_cast<int>(triangles.size());
    const int nr = static_cast<int>(retained.size());

    // ---- LST nodes: the vertices, then one midside node per unique edge
    std::map<std::pair<int, int>, int> edgeNode;
    std::vector<std::array<int, 6>> conn(nt);
    int nn = nv;
    for (int e = 0; e < nt; ++e) {
        const auto& tri = triangles[e];
        conn[e][0] = tri[0]; conn[e][1] = tri[1]; conn[e][2] = tri[2];
        const int ends[3][2] = {{tri[0], tri[1]}, {tri[1], tri[2]}, {tri[2], tri[0]}};
        for (int k = 0; k < 3; ++k) {
            std::pair<int, int> key(std::min(ends[k][0], ends[k][1]), std::max(ends[k][0], ends[k][1]));
            auto it = edgeNode.find(key);
            if (it == edgeNode.end())
                it = edgeNode.emplace(key, nn++).first;
            conn[e][3 + k] = it->second;
        }
    }

    // ---- DOF partition: retained first (slot order), then every other node
    // that belongs to a triangle; nodes in no triangle are dropped.
    std::vector<char> used(nn, 0);
    for (const auto& cn : conn)
        for (int a : cn) used[a] = 1;
    std::vector<int> nodeSlot(nn, -1);        // node -> position in the reduced list
    Condensed out;
    out.nr = nr;
    out.nTri = nt;
    for (int r = 0; r < nr; ++r) {
        int v = retained[r];
        if (v < 0 || v >= nv || !used[v]) {
            out.isolated.push_back(r);
            continue;
        }
        if (nodeSlot[v] >= 0)
            throw std::invalid_argument("two retained nodes share one mesh vertex");
        nodeSlot[v] = r;
    }
    int ni = 0;
    std::vector<int> interiorOf(nn, -1);
    for (int a = 0; a < nn; ++a)
        if (used[a] && nodeSlot[a] < 0)
            interiorOf[a] = ni++;
    const int nrd = 2 * nr, nid = 2 * ni;
    // global reduced dof of node a, component k: retained [0, nrd), interior [nrd, ...)
    auto dofOf = [&](int a, int k) -> int {
        return (nodeSlot[a] >= 0) ? 2 * nodeSlot[a] + k : nrd + 2 * interiorOf[a] + k;
    };

    // ---- assemble K and M (LST, exact quadrature)
    std::vector<Trip> kt, mt;
    kt.reserve(static_cast<size_t>(nt) * 144);
    if (rhoT > 0.0) mt.reserve(static_cast<size_t>(nt) * 72);
    out.triArea.assign(nt, 0.0);
    // midpoint rule (exact for the quadratic B'DB)
    static const double QK[3][3] = {{0.5, 0.5, 0.0}, {0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}};
    std::vector<std::array<double, 3>> bAll(nt), cAll(nt);
    std::vector<double> A2All(nt);
    for (int e = 0; e < nt; ++e) {
        const Point2D& p0 = vertices[conn[e][0]];
        const Point2D& p1 = vertices[conn[e][1]];
        const Point2D& p2 = vertices[conn[e][2]];
        double b[3] = {p1.y - p2.y, p2.y - p0.y, p0.y - p1.y};
        double c[3] = {p2.x - p1.x, p0.x - p2.x, p1.x - p0.x};
        double A2 = (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y);
        if (A2 <= 0.0)
            throw std::invalid_argument("degenerate or clockwise triangle in the diaphragm mesh");
        bAll[e] = {b[0], b[1], b[2]};
        cAll[e] = {c[0], c[1], c[2]};
        A2All[e] = A2;
        double area = 0.5 * A2;
        out.triArea[e] = area;

        double ke[12][12] = {};
        for (const auto& q : QK) {
            double B[3][12];
            detail::lst_B(q, b, c, A2, B);
            double DB[3][12];
            for (int r = 0; r < 3; ++r)
                for (int j = 0; j < 12; ++j)
                    DB[r][j] = D[3 * r] * B[0][j] + D[3 * r + 1] * B[1][j] + D[3 * r + 2] * B[2][j];
            double w = t * area / 3.0;
            for (int i = 0; i < 12; ++i)
                for (int j = 0; j < 12; ++j)
                    ke[i][j] += w * (B[0][i] * DB[0][j] + B[1][i] * DB[1][j] + B[2][i] * DB[2][j]);
        }
        int gd[12];
        for (int a = 0; a < 6; ++a) {
            gd[2 * a] = dofOf(conn[e][a], 0);
            gd[2 * a + 1] = dofOf(conn[e][a], 1);
        }
        for (int i = 0; i < 12; ++i)
            for (int j = 0; j < 12; ++j)
                if (ke[i][j] != 0.0)
                    kt.emplace_back(gd[i], gd[j], ke[i][j]);

        if (rhoT > 0.0) {
            double me[6][6] = {};
            for (std::size_t g = 0; g < LadrunoMembraneKernel::TRI7_W.size(); ++g) {
                const auto& L = LadrunoMembraneKernel::TRI7_BARY[g];
                double Lq[3] = {L[0], L[1], L[2]};
                double N[6];
                detail::lst_N(Lq, N);
                double w = rhoT * area * LadrunoMembraneKernel::TRI7_W[g];
                for (int i = 0; i < 6; ++i)
                    for (int j = 0; j < 6; ++j)
                        me[i][j] += w * N[i] * N[j];
            }
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j)
                    for (int k = 0; k < 2; ++k)
                        mt.emplace_back(gd[2 * i + k], gd[2 * j + k], me[i][j]);
        }
    }

    const int ndof = nrd + nid;
    Eigen::SparseMatrix<double> K(ndof, ndof), M(ndof, ndof);
    K.setFromTriplets(kt.begin(), kt.end());
    if (rhoT > 0.0)
        M.setFromTriplets(mt.begin(), mt.end());

    // ---- condensation: X = -K_ii^-1 K_ir  (nid x nrd)
    Eigen::MatrixXd X(nid, nrd);
    Eigen::SparseMatrix<double> Kii = K.block(nrd, nrd, nid, nid);
    Eigen::SparseMatrix<double> Kir = K.block(nrd, 0, nid, nrd);
    if (nid > 0) {
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
        ldlt.compute(Kii);
        if (ldlt.info() != Eigen::Success)
            throw std::runtime_error("diaphragm condensation: K_ii factorization failed");
        Eigen::MatrixXd rhs = -Eigen::MatrixXd(Kir);
        X = ldlt.solve(rhs);
        if (ldlt.info() != Eigen::Success)
            throw std::runtime_error("diaphragm condensation: K_ii solve failed");
    }

    Eigen::MatrixXd Kc = Eigen::MatrixXd(K.block(0, 0, nrd, nrd));
    if (nid > 0)
        Kc += Eigen::MatrixXd(K.block(0, nrd, nrd, nid)) * X;
    Kc = 0.5 * (Kc + Kc.transpose()).eval();
    out.K.assign(Kc.data(), Kc.data() + static_cast<size_t>(nrd) * nrd); // symmetric: col = row major

    out.M.assign(static_cast<size_t>(nrd) * nrd, 0.0);
    if (rhoT > 0.0) {
        Eigen::MatrixXd Mrr = Eigen::MatrixXd(M.block(0, 0, nrd, nrd));
        if (nid > 0) {
            Eigen::MatrixXd MirX = Eigen::MatrixXd(M.block(0, nrd, nrd, nid)) * X;
            Eigen::MatrixXd MiiX = M.block(nrd, nrd, nid, nid) * X;
            Mrr += MirX + MirX.transpose() + X.transpose() * MiiX;
        }
        Mrr = 0.5 * (Mrr + Mrr.transpose()).eval();
        out.M.assign(Mrr.data(), Mrr.data() + static_cast<size_t>(nrd) * nrd);
    }

    // ---- strain recovery at triangle centroids: eps = B(centroid) R_e u_r
    out.strainOp.assign(static_cast<size_t>(nt) * 3 * nrd, 0.0);
    out.meanStrainOp.assign(static_cast<size_t>(3) * nrd, 0.0);
    double totalArea = 0.0;
    const double Lc[3] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    for (int e = 0; e < nt; ++e) {
        double B[3][12];
        detail::lst_B(Lc, bAll[e].data(), cAll[e].data(), A2All[e], B);
        double* op = &out.strainOp[static_cast<size_t>(e) * 3 * nrd];
        for (int a = 0; a < 6; ++a) {
            int node = conn[e][a];
            for (int k = 0; k < 2; ++k) {
                int col = 2 * a + k;
                if (nodeSlot[node] >= 0) {
                    int rd = 2 * nodeSlot[node] + k;
                    for (int r = 0; r < 3; ++r)
                        op[r * nrd + rd] += B[r][col];
                } else {
                    int id = 2 * interiorOf[node] + k;
                    for (int r = 0; r < 3; ++r) {
                        double brc = B[r][col];
                        if (brc == 0.0) continue;
                        for (int j = 0; j < nrd; ++j)
                            op[r * nrd + j] += brc * X(id, j);
                    }
                }
            }
        }
        double area = out.triArea[e];
        totalArea += area;
        for (int i = 0; i < 3 * nrd; ++i)
            out.meanStrainOp[i] += area * op[i];
    }
    if (totalArea > 0.0)
        for (double& v : out.meanStrainOp) v /= totalArea;

    return out;
}

} // namespace LadrunoMembraneDiaphragm

#endif // LadrunoMembraneDiaphragm_h
