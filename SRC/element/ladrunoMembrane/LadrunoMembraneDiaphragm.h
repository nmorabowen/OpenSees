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

#include <algorithm>
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

// WP-208: a connection line between two consecutive retained nodes (a loop
// segment of the outer polygon or a hole, split at every retained node on it,
// or a supporting frame whose ends are element nodes). Every mesh node strictly
// inside it (vertex or LST midside) is SLAVED to the linear interpolation of
// the two ends: the slab edge between connection nodes is straight, as in the
// linear-edge element it is connected to (a wall or beam mesh, or the ETABS
// slab mesh whose boundary nodes the loop vertices are). Without it the edge
// bowed between connection nodes and the panel failed the constant-strain
// patch test (reaction 0.74 of exact, worse with refinement).
struct LinkLine {
    Point2D a, b;
    int slotA = -1, slotB = -1;
};

// vertices/triangles: the panel mesh (ccw triangles). retained[r] is the mesh
// vertex of retained slot r (-1: none — the slot gets no stiffness). D is the
// 3x3 plane-stress modulus (row-major, modifiers already applied), t the
// thickness, rhoT the mass per unit area (0: no mass), links the connection
// lines (WP-208, see LinkLine).
inline Condensed condense_lst(const std::vector<Point2D>& vertices,
                              const std::vector<std::array<int, 3>>& triangles,
                              const std::vector<int>& retained,
                              const double D[9], double t, double rhoT,
                              const std::vector<LinkLine>& links = {})
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
    // node positions (midsides at the edge midpoints)
    std::vector<Point2D> pos(nn);
    for (int v = 0; v < nv; ++v) pos[v] = vertices[v];
    for (const auto& kv : edgeNode) {
        const Point2D& pa = vertices[kv.first.first];
        const Point2D& pb = vertices[kv.first.second];
        pos[kv.second] = {0.5 * (pa.x + pb.x), 0.5 * (pa.y + pb.y)};
    }
    // WP-208: slave non-retained nodes lying strictly inside a connection line
    std::vector<int> slaveLink(nn, -1);
    std::vector<double> slaveS(nn, 0.0);
    double scale = 0.0;
    for (const auto& p : vertices) scale = std::max({scale, std::abs(p.x), std::abs(p.y)});
    const double tolD = 1e-9 * std::max(scale, 1.0);
    for (int a = 0; a < nn; ++a) {
        if (!used[a] || nodeSlot[a] >= 0) continue;
        for (std::size_t l = 0; l < links.size(); ++l) {
            const LinkLine& L = links[l];
            if (L.slotA < 0 || L.slotB < 0 || L.slotA >= nr || L.slotB >= nr) continue;
            if (retained[L.slotA] < 0 || retained[L.slotB] < 0) continue;
            double dx = L.b.x - L.a.x, dy = L.b.y - L.a.y;
            double len2 = dx * dx + dy * dy;
            if (len2 <= 0.0) continue;
            double s = ((pos[a].x - L.a.x) * dx + (pos[a].y - L.a.y) * dy) / len2;
            if (s <= 1e-9 || s >= 1.0 - 1e-9) continue;
            double ex = pos[a].x - (L.a.x + s * dx), ey = pos[a].y - (L.a.y + s * dy);
            if (std::hypot(ex, ey) > tolD) continue;
            slaveLink[a] = static_cast<int>(l);
            slaveS[a] = s;
            break;
        }
    }
    int ni = 0;
    std::vector<int> interiorOf(nn, -1);
    for (int a = 0; a < nn; ++a)
        if (used[a] && nodeSlot[a] < 0 && slaveLink[a] < 0)
            interiorOf[a] = ni++;
    const int nrd = 2 * nr, nid = 2 * ni;
    // node a, component k as a combination of reduced dofs: retained [0, nrd),
    // interior [nrd, ...); a slaved node is (1-s) end A + s end B (retained)
    struct Term { int dof; double w; };
    auto termsOf = [&](int a, int k, Term out2[2]) -> int {
        if (nodeSlot[a] >= 0) { out2[0] = {2 * nodeSlot[a] + k, 1.0}; return 1; }
        if (slaveLink[a] >= 0) {
            const LinkLine& L = links[slaveLink[a]];
            out2[0] = {2 * L.slotA + k, 1.0 - slaveS[a]};
            out2[1] = {2 * L.slotB + k, slaveS[a]};
            return 2;
        }
        out2[0] = {nrd + 2 * interiorOf[a] + k, 1.0};
        return 1;
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
        Term gt[12][2];
        int gn[12];
        for (int a = 0; a < 6; ++a)
            for (int k = 0; k < 2; ++k)
                gn[2 * a + k] = termsOf(conn[e][a], k, gt[2 * a + k]);
        for (int i = 0; i < 12; ++i)
            for (int j = 0; j < 12; ++j) {
                if (ke[i][j] == 0.0) continue;
                for (int p = 0; p < gn[i]; ++p)
                    for (int q = 0; q < gn[j]; ++q)
                        kt.emplace_back(gt[i][p].dof, gt[j][q].dof, ke[i][j] * gt[i][p].w * gt[j][q].w);
            }

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
                    for (int k = 0; k < 2; ++k) {
                        int ii = 2 * i + k, jj = 2 * j + k;
                        for (int p = 0; p < gn[ii]; ++p)
                            for (int q = 0; q < gn[jj]; ++q)
                                mt.emplace_back(gt[ii][p].dof, gt[jj][q].dof,
                                                me[i][j] * gt[ii][p].w * gt[jj][q].w);
                    }
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
                Term tm[2];
                int nterm = termsOf(node, k, tm);
                for (int p = 0; p < nterm; ++p) {
                    if (tm[p].dof < nrd) {
                        for (int r = 0; r < 3; ++r)
                            op[r * nrd + tm[p].dof] += B[r][col] * tm[p].w;
                    } else {
                        int id = tm[p].dof - nrd;
                        for (int r = 0; r < 3; ++r) {
                            double brc = B[r][col] * tm[p].w;
                            if (brc == 0.0) continue;
                            for (int j = 0; j < nrd; ++j)
                                op[r * nrd + j] += brc * X(id, j);
                        }
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
