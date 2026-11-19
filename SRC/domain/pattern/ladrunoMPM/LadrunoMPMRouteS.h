/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno: WP-232l route S -- the semi-implicit u-p step of mode F (tests/mpm_twophase_interface.md
// §1 "Route S (mode F) step", §5 P11, §7a Q4, §7b P5a and P4 amended; wp-briefs "232l" rulings
// 9064acca3 (1, 2), 3d3699a50, c84f71fb3, bb1421559). Included by LadrunoMPMCore.h only.
//
// One step, explicit skeleton and implicit pressure (IMEX):
//   1. P2G with the per-component folded supports (upSupport): lumped W mass, momentum, the
//      effective-stress force, gravity on the mixture, tractions; plates as one normal unknown U.
//   2. Predictor v* = v^ + dt (f' + f_ext)/m (local damping as mode F), then the velocity
//      conditions; M^-1 = 1/m on a free component, 0 on a fixed, prescribed, stray or inactive one.
//   3. Pressure on the p-active, non-drained nodes:
//        [S* + dt H + dt^2 Q^T M^-1 Q] p = S* p^ - dt Q^T v* + dt f_seep,
//      drained p0 moved to the right side; Jacobi-PCG rel 1e-13, maxit 4 n + 1000, warm start p^.
//      p^ is P5a's consistent L2 projection over the p-active nodes (probe-guarded PCG).
//   4. Correction v1 = v* + dt M^-1 Q p on free components; a = (v1 - v^)/dt.
//   5. MUSL: particle v, x from v1; second projection; dF = I + dt L (Delta eps = sym(dt L) under
//      `small`); sigma'.
//   6. G2P of p (P5a: FLIP on p - p^, blend alpha_p); commit.

#ifndef LadrunoMPMRouteS_h
#define LadrunoMPMRouteS_h

#include <map>
#include <set>

namespace ladruno_mpm {

// Jacobi-PCG on a symmetric matrix stored as row maps. relTol is on ||r|| / ||b||; guard: refuse a
// curvature d^T A d <= 1e-14 max diag(A) d^T d (P5a's frozen guard). x holds the warm start.
// Returns the iteration count (>= 0), -1 when not converged within maxit, -2 on the guard.
inline int routeSPcg(const std::vector<std::map<int, double>> &A, const std::vector<double> &b,
                     std::vector<double> &x, double relTol, bool guard)
{
  const int n = (int)A.size();
  std::vector<double> diag(n, 0.0), r(n), z(n), d(n), Ad(n);
  double dmax = 0.0;
  for (int i = 0; i < n; i++) {
    const auto it = A[i].find(i);
    diag[i] = (it != A[i].end()) ? it->second : 0.0;
    if (!(diag[i] > 0.0)) return -2;
    dmax = std::max(dmax, diag[i]);
  }
  if ((int)x.size() != n) x.assign(n, 0.0);
  double bn = 0.0;
  for (int i = 0; i < n; i++) bn += b[i] * b[i];
  bn = std::sqrt(bn);
  for (int i = 0; i < n; i++) {
    double s = 0.0;
    for (const auto &kv : A[i]) s += kv.second * x[kv.first];
    r[i] = b[i] - s;
  }
  auto norm = [&](const std::vector<double> &v) {
    double s = 0.0;
    for (double q : v) s += q * q;
    return std::sqrt(s);
  };
  if (!(bn > 0.0)) {
    if (!(norm(r) > 0.0)) return 0;
    bn = 1.0;
  }
  if (norm(r) <= relTol * bn) return 0;
  double rz = 0.0;
  for (int i = 0; i < n; i++) { z[i] = r[i] / diag[i]; d[i] = z[i]; rz += r[i] * z[i]; }
  const int maxit = 4 * n + 1000;
  for (int it = 0; it < maxit; it++) {
    double dAd = 0.0, dd = 0.0;
    for (int i = 0; i < n; i++) {
      double s = 0.0;
      for (const auto &kv : A[i]) s += kv.second * d[kv.first];
      Ad[i] = s;
      dAd += d[i] * s;
      dd += d[i] * d[i];
    }
    if (guard && !(dAd > 1e-14 * dmax * dd)) return -2;
    if (!(dAd > 0.0)) return -2;
    const double al = rz / dAd;
    for (int i = 0; i < n; i++) { x[i] += al * d[i]; r[i] -= al * Ad[i]; }
    if (norm(r) <= relTol * bn) return it + 1;
    double rz1 = 0.0;
    for (int i = 0; i < n; i++) { z[i] = r[i] / diag[i]; rz1 += r[i] * z[i]; }
    const double be = rz1 / rz;
    rz = rz1;
    for (int i = 0; i < n; i++) d[i] = z[i] + be * d[i];
  }
  return -1;
}

inline int Engine::advanceUP(double dt, double lambda, MaterialPoints &mats, UPStep &u, std::string &err)
{
  const int N = P.size(), nd = grid.ndm, nvv = nv(), nc = nd + 1, nn = nd * nd;
  if (N == 0) return 0;
  if (!bodyForce.empty() && (long)bodyForce.size() != (long)nd * N) {
    err = "the particle body force holds " + std::to_string((long long)bodyForce.size()) +
          " values for " + std::to_string(N) + " particles (ndm per particle): set it again;"
          " step refused, nothing changed";
    return -1;
  }
  if ((int)u.tp.size() != N || (int)u.alphaSt.size() != N) {
    err = "internal: the route-S constants do not match the particles; step refused";
    return -1;
  }
  bool up = false;
  for (int p = 0; p < N; p++) up = up || u.tp[p] != 0;

  // ---- 1. supports (P9 per-component folds, P10 plate entries) and P2G ---------------------
  std::vector<long> sS(1, 0), sNode;
  std::vector<double> sW, sdW, sSc, sdSc;
  for (int p = 0; p < N; p++) {
    upSupport(grid, par.basis, fc, &P.x[nd * p], &P.lp[nd * p], sNode, sW, sdW, sSc, sdSc, 0);
    sS.push_back((long)sNode.size());
  }
  std::unordered_map<long, int> loc;          // lattice node -> local index
  std::vector<long> nodes;
  for (long x : sNode)
    if (x >= 0 && loc.find(x) == loc.end()) { loc[x] = (int)nodes.size(); nodes.push_back(x); }
  const int nl = (int)nodes.size();
  // per entry: the local index (>= 0), or -(1 + f) for plate face f
  std::vector<int> eIdx(sNode.size());
  for (size_t k = 0; k < sNode.size(); k++)
    eIdx[k] = (sNode[k] >= 0) ? loc[sNode[k]] : -(1 + (int)(-sNode[k] - 2));
  std::vector<double> m(nl, 0.0), tpm(nl, 0.0), gp(nd * nl, 0.0), gf(nd * nl, 0.0);
  double mU[6] = {0, 0, 0, 0, 0, 0}, pU[6] = {0, 0, 0, 0, 0, 0}, fU[6] = {0, 0, 0, 0, 0, 0};
  double b[3];
  for (int d = 0; d < 3; d++) b[d] = lambda * par.gravity[d];
  for (int p = 0; p < N; p++) {
    const double mp = P.m[p], Vp = P.V[p];
    const double *bfp = bodyForce.empty() ? nullptr : &bodyForce[nd * p];
    double sg[3][3];
    voigtToTensor(&P.sig[nvv * p], sg);
    for (long k = sS[p]; k < sS[p + 1]; k++) {
      const int i = eIdx[k];
      if (i >= 0) {
        m[i] += sW[k] * mp;
        if (u.tp[p] != 0) tpm[i] += sW[k] * mp;
      } else {
        mU[-i - 1] += sW[k] * mp;
      }
      for (int c = 0; c < nd; c++) {
        const double S = sSc[nc * k + c];
        const double *g = &sdSc[(nc * k + c) * nd];
        double s = 0.0;
        for (int j = 0; j < nd; j++) s += sg[c][j] * g[j];
        const double mom = S * mp * P.v[nd * p + c];
        const double frc = S * mp * (bfp ? b[c] + bfp[c] : b[c]) - Vp * s;
        if (i >= 0) { gp[nd * i + c] += mom; gf[nd * i + c] += frc; }
        else if (c == (-i - 1) / 2) { pU[-i - 1] += mom; fU[-i - 1] += frc; }
      }
    }
  }
  for (const UPStep::NodeForce &nf : u.nodeForce) {
    const auto it = loc.find(nf.node);
    if (it == loc.end()) continue;
    for (int c = 0; c < nd; c++) gf[nd * it->second + c] += nf.f[c];
  }
  for (int f = 0; f < 2 * nd; f++) fU[f] += u.plateLoad[f];

  // ---- 2. predictor and the velocity conditions; M^-1 (§7a Q4) -----------------------------
  const double mmax = massRef();
  const double mcut = par.masscut * mmax;
  std::vector<char> act(nl, 0);
  std::vector<double> vhat(nd * nl, 0.0), vst(nd * nl, 0.0), minv(nd * nl, 0.0);
  const double SENT = 1.2345678e300;
  for (int i = 0; i < nl; i++) {
    act[i] = (m[i] > 0.0 && imageCutoffMass(grid, par.basis, fc, nodes[i], m[i]) > mcut) ? 1 : 0;
    if (!act[i]) continue;
    double v[3], probe[3] = {SENT, SENT, SENT};
    for (int c = 0; c < nd; c++) {
      vhat[nd * i + c] = gp[nd * i + c] / m[i];
      double f = gf[nd * i + c];
      const double vn = vhat[nd * i + c];
      const double sgn = (vn > 0.0) ? 1.0 : ((vn < 0.0) ? -1.0 : 0.0);
      f -= par.damping * std::fabs(f) * sgn;
      v[c] = vn + dt * f / m[i];
    }
    applyConditions(nodes[i], v, true);
    applyConditions(nodes[i], probe, true);
    int ijk[3];
    grid.extIJK(nodes[i], ijk);
    for (int c = 0; c < nd; c++) {
      bool stray = false;                  // P10: a plate face row's normal component
      for (int f = 0; f < 2 * nd; f++)
        if (fc.plateOn[f] && c == f / 2 && ijk[c] == ((f % 2) ? grid.n[c] : 0)) stray = true;
      if (stray) v[c] = 0.0;
      vst[nd * i + c] = v[c];
      minv[nd * i + c] = (stray || probe[c] != SENT) ? 0.0 : 1.0 / m[i];
    }
  }
  double vhatU[6] = {0, 0, 0, 0, 0, 0}, vstU[6] = {0, 0, 0, 0, 0, 0}, minvU[6] = {0, 0, 0, 0, 0, 0};
  for (int f = 0; f < 2 * nd; f++) {
    if (!fc.plateOn[f] || !(mU[f] > 0.0)) continue;
    vhatU[f] = pU[f] / mU[f];
    double F = fU[f];
    const double sgn = (vhatU[f] > 0.0) ? 1.0 : ((vhatU[f] < 0.0) ? -1.0 : 0.0);
    F -= par.damping * std::fabs(F) * sgn;
    vstU[f] = vhatU[f] + dt * F / mU[f];
    minvU[f] = 1.0 / mU[f];
  }
  std::vector<double> v1(vst);
  double v1U[6];
  for (int f = 0; f < 6; f++) v1U[f] = vstU[f];

  // ---- 3. the pressure (two-phase epochs) ----------------------------------------------------
  std::vector<double> SpF, dSpF;
  std::vector<int> pidx(nl, -1);
  std::vector<int> pnode;                    // p unknown -> local node
  std::vector<double> phat, pa;              // per p unknown: projection, solution
  std::vector<double> tpw(P.pw);
  if (up) {
    // p-activity (3d3699a50, c84f71fb3; P4 amended): two-phase nodal mass > pCut m_max, and
    // some two-phase particle with S^p > 0 there
    std::vector<char> reach(nl, 0);
    for (int p = 0; p < N; p++) {
      if (u.tp[p] == 0) continue;
      for (long k = sS[p]; k < sS[p + 1]; k++)
        if (eIdx[k] >= 0 && sSc[nc * k + nd] > 0.0) reach[eIdx[k]] = 1;
    }
    for (int i = 0; i < nl; i++)
      if (reach[i] && imageCutoffMass(grid, par.basis, fc, nodes[i], tpm[i]) > u.pCut * mmax) {
        pidx[i] = (int)pnode.size();
        pnode.push_back(i);
      }
    const int np = (int)pnode.size();
    // S^p renormalised over the p-active nodes on the RAW weights; a folded-away entry (S^p = 0)
    // counts as retained (bb1421559)
    SpF.assign(sNode.size(), 0.0);
    dSpF.assign(nd * sNode.size(), 0.0);
    for (int p = 0; p < N; p++) {
      if (u.tp[p] == 0) continue;
      bool allAct = true, anyAct = false;
      double sum = 0.0, gsum[3] = {0.0, 0.0, 0.0};
      for (long k = sS[p]; k < sS[p + 1]; k++) {
        if (eIdx[k] < 0) continue;
        const bool pa = pidx[eIdx[k]] >= 0;
        if (pa) anyAct = true;
        if (pa || !(sSc[nc * k + nd] > 0.0)) {
          sum += sW[k];
          for (int a = 0; a < nd; a++) gsum[a] += sdW[nd * k + a];
        } else {
          allAct = false;
        }
      }
      if (!anyAct || !(sum > 0.0)) {
        err = "two-phase particle " + std::to_string((long long)P.id[p]) + " has no p-active node (its"
              " two-phase nodal mass is below `mpm pressureCutoff` everywhere); step refused, nothing"
              " changed (interface 7a Q4)";
        return -1;
      }
      for (long k = sS[p]; k < sS[p + 1]; k++) {
        if (eIdx[k] < 0 || pidx[eIdx[k]] < 0) continue;
        double sp = sSc[nc * k + nd], gg[3];
        for (int a = 0; a < nd; a++) gg[a] = sdSc[(nc * k + nd) * nd + a];
        if (!allAct) {
          sp /= sum;
          for (int a = 0; a < nd; a++) gg[a] = (gg[a] - sp * gsum[a]) / sum;
        }
        SpF[k] = sp;
        for (int a = 0; a < nd; a++) dSpF[nd * k + a] = gg[a];
      }
    }
    // drained nodes (P8): a p-active node on a `pfix` face or its ghost layer
    std::vector<char> drained(np, 0);
    std::vector<double> p0(np, 0.0);
    for (int a = 0; a < np; a++) {
      int ijk[3];
      grid.extIJK(nodes[pnode[a]], ijk);
      for (int f = 0; f < 2 * nd; f++) {
        const int d = f / 2;
        const bool on = (f % 2 == 0) ? ijk[d] <= 0 : ijk[d] >= grid.n[d];
        if (!on || !fc.pfixOn[f]) continue;
        if (drained[a] && p0[a] != fc.pfixVal[f]) {
          err = "two `pfix` faces prescribe different p on one node (a corner): refused";
          return -1;
        }
        drained[a] = 1;
        p0[a] = fc.pfixVal[f];
      }
    }
    // P5a: the consistent L2 projection M p^ = P^T V p over the p-active nodes (pinned included)
    std::vector<std::map<int, double>> Mp(np);
    std::vector<double> rp(np, 0.0);
    for (int p = 0; p < N; p++) {
      if (u.tp[p] == 0) continue;
      const double V = P.V[p];
      for (long k = sS[p]; k < sS[p + 1]; k++) {
        if (SpF[k] == 0.0) continue;
        const int a = pidx[eIdx[k]];
        rp[a] += V * SpF[k] * P.pw[p];
        for (long l = sS[p]; l < sS[p + 1]; l++)
          if (SpF[l] != 0.0) Mp[a][pidx[eIdx[l]]] += V * SpF[k] * SpF[l];
      }
    }
    const std::string singular =
        "the consistent pressure projection (interface 7b, P5a) is singular or did not reach rel 1e-14:"
        " the particles no longer span the p-active nodes' shape functions (e.g. one particle per cell"
        " on a B-spline grid); step refused, nothing changed";
    phat.assign(np, 0.0);
    if (np > 0) {
      if (routeSPcg(Mp, rp, phat, 1e-14, true) < 0) { err = singular; return -1; }
      std::vector<double> bq(np), xq;                     // the probe (guard (a), ruled)
      for (int a = 0; a < np; a++) bq[a] = std::sin(1.7 * a + 0.3);
      if (routeSPcg(Mp, bq, xq, 1e-14, true) < 0) { err = singular + " (the probe solve)"; return -1; }
    }
    // the blocks (kernel helpers' forms, §1): S* = S + alpha_st H~, H, f_seep, Q
    std::vector<std::map<int, double>> A(np), Ss(np);
    std::vector<double> rhs(np, 0.0);
    std::map<int, std::map<int, double>> Q;          // u row: nd * local + c, or -(1 + f) plate
    for (int p = 0; p < N; p++) {
      const TwoPhaseProps *tp = u.tp[p];
      if (tp == 0) continue;
      const double V = P.V[p], iQ = tp->oneOverQbar(), ast = u.alphaSt[p];
      double drv[3] = {0.0, 0.0, 0.0};
      for (int d = 0; d < nd; d++) drv[d] = lambda * (tp->fluidBodyGiven ? tp->fluidBody[d] : par.gravity[d]);
      for (long k = sS[p]; k < sS[p + 1]; k++) {
        if (SpF[k] == 0.0) continue;
        const int a = pidx[eIdx[k]];
        const double *ga = &dSpF[nd * k];
        double fs = 0.0;
        for (int d = 0; d < nd; d++) fs += ga[d] * tp->kbar[d] * tp->rhoF * drv[d];
        rhs[a] += dt * V * fs;
        for (long l = sS[p]; l < sS[p + 1]; l++) {
          if (SpF[l] == 0.0) continue;
          const int c = pidx[eIdx[l]];
          const double *gc = &dSpF[nd * l];
          double gg = 0.0, hk = 0.0;
          for (int d = 0; d < nd; d++) { gg += ga[d] * gc[d]; hk += ga[d] * tp->kbar[d] * gc[d]; }
          const double sst = V * (iQ * SpF[k] * SpF[l] + ast * gg);
          Ss[a][c] += sst;
          A[a][c] += sst + dt * V * hk;
        }
        // Q_(b,i),a = alpha V d_i S'_{i,b} S^p_a over every u entry of the particle
        for (long l = sS[p]; l < sS[p + 1]; l++)
          for (int i = 0; i < nd; i++) {
            const int e = eIdx[l];
            if (e < 0 && i != (-e - 1) / 2) continue;
            const double g = sdSc[(nc * l + i) * nd + i];
            if (g == 0.0) continue;
            const int row = (e >= 0) ? nd * e + i : e;
            Q[row][a] += tp->alpha * V * g * SpF[k];
          }
      }
    }
    // + dt^2 Q^T M^-1 Q;  rhs += S* p^ - dt Q^T v*
    for (const auto &rw : Q) {
      const int row = rw.first;
      const double mi = (row >= 0) ? minv[row] : minvU[-row - 1];
      const double vs = (row >= 0) ? vst[row] : vstU[-row - 1];
      for (const auto &qa : rw.second) {
        rhs[qa.first] -= dt * qa.second * vs;
        if (mi != 0.0)
          for (const auto &qc : rw.second) A[qa.first][qc.first] += dt * dt * mi * qa.second * qc.second;
      }
    }
    for (int a = 0; a < np; a++)
      for (const auto &kv : Ss[a]) rhs[a] += kv.second * phat[kv.first];
    // the free system: drained p0 to the right side
    std::vector<int> fidx(np, -1);
    std::vector<int> fre;
    for (int a = 0; a < np; a++) if (!drained[a]) { fidx[a] = (int)fre.size(); fre.push_back(a); }
    const int nf = (int)fre.size();
    std::vector<std::map<int, double>> Af(nf);
    std::vector<double> bf(nf), xf(nf);
    for (int r = 0; r < nf; r++) {
      const int a = fre[r];
      bf[r] = rhs[a];
      xf[r] = phat[a];                                    // warm start p^
      for (const auto &kv : A[a]) {
        if (drained[kv.first]) bf[r] -= kv.second * p0[kv.first];
        else Af[r][fidx[kv.first]] += kv.second;
      }
    }
    pa.assign(np, 0.0);
    if (nf > 0) {
      const int it = routeSPcg(Af, bf, xf, 1e-13, false);
      if (it < 0) {
        err = "the route-S pressure solve (interface §1 step 3) did not reach rel 1e-13 in 4 n + 1000"
              " Jacobi-PCG iterations; step refused, nothing changed";
        return -1;
      }
      u.pcgIter = it;
    }
    for (int a = 0; a < np; a++) pa[a] = drained[a] ? p0[a] : xf[fidx[a]];
    // ---- 4. correction v1 = v* + dt M^-1 Q p ---------------------------------------------
    for (const auto &rw : Q) {
      const int row = rw.first;
      const double mi = (row >= 0) ? minv[row] : minvU[-row - 1];
      if (mi == 0.0) continue;
      double s = 0.0;
      for (const auto &qa : rw.second) s += qa.second * pa[qa.first];
      if (row >= 0) v1[row] += dt * mi * s;
      else v1U[-row - 1] += dt * mi * s;
    }
    // ---- 6 (trial). G2P of p: P5a, FLIP on p - p^ with the blend alpha_p --------------------
    for (int p = 0; p < N; p++) {
      if (u.tp[p] == 0) continue;
      double inc = 0.0, pic = 0.0;
      for (long k = sS[p]; k < sS[p + 1]; k++) {
        if (SpF[k] == 0.0) continue;
        const int a = pidx[eIdx[k]];
        inc += SpF[k] * (pa[a] - phat[a]);
        pic += SpF[k] * pa[a];
      }
      tpw[p] = u.pBlend * (P.pw[p] + inc) + (1.0 - u.pBlend) * pic;
    }
  }
  std::vector<double> ag(nd * nl, 0.0);
  double aU[6] = {0, 0, 0, 0, 0, 0};
  for (int i = 0; i < nl; i++)
    if (act[i])
      for (int c = 0; c < nd; c++) ag[nd * i + c] = (v1[nd * i + c] - vhat[nd * i + c]) / dt;
  for (int f = 0; f < 6; f++) if (minvU[f] != 0.0) aU[f] = (v1U[f] - vhatU[f]) / dt;

  // ---- 5. MUSL particle update (trial) ---------------------------------------------------------
  tx.assign(P.x.begin(), P.x.end());
  tv.assign(P.v.begin(), P.v.end());
  for (int p = 0; p < N; p++) {
    double sa[3] = {0, 0, 0}, sv[3] = {0, 0, 0};
    for (long k = sS[p]; k < sS[p + 1]; k++) {
      const int i = eIdx[k];
      for (int c = 0; c < nd; c++) {
        const double S = sSc[nc * k + c];
        if (S == 0.0) continue;
        if (i >= 0) { sa[c] += S * ag[nd * i + c]; sv[c] += S * v1[nd * i + c]; }
        else if (c == (-i - 1) / 2) { sa[c] += S * aU[-i - 1]; sv[c] += S * v1U[-i - 1]; }
      }
    }
    for (int c = 0; c < nd; c++) {
      tv[nd * p + c] = par.flip * (P.v[nd * p + c] + dt * sa[c]) + (1.0 - par.flip) * sv[c];
      tx[nd * p + c] = P.x[nd * p + c] + dt * sv[c];
    }
  }
  {
    int worst = -1;
    for (int p = 0; p < N; p++)
      for (int a = 0; a < nd; a++)
        if (!admissibleAxis(grid, par.basis, a, tx[nd * p + a], P.lp[nd * p + a])) {
          if (worst < 0 || P.id[p] < P.id[worst]) worst = p;
          break;
        }
    if (worst >= 0) {
      err = "the step would carry particle " + std::to_string((long long)P.id[worst]) +
            " out of the grid (interface A7); step refused, nothing changed";
      return -1;
    }
  }
  // second projection, conditioned (no walls: refused with route S)
  std::vector<double> gvs(nd * nl, 0.0);
  double gvsU[6] = {0, 0, 0, 0, 0, 0};
  for (int p = 0; p < N; p++)
    for (long k = sS[p]; k < sS[p + 1]; k++) {
      const int i = eIdx[k];
      for (int c = 0; c < nd; c++) {
        const double w = sSc[nc * k + c] * P.m[p] * tv[nd * p + c];
        if (i >= 0) gvs[nd * i + c] += w;
        else if (c == (-i - 1) / 2) gvsU[-i - 1] += w;
      }
    }
  for (int i = 0; i < nl; i++) {
    double vs[3] = {0, 0, 0};
    if (act[i]) {
      for (int c = 0; c < nd; c++) vs[c] = gvs[nd * i + c] / m[i];
      applyConditions(nodes[i], vs, false);
      int ijk[3];
      grid.extIJK(nodes[i], ijk);
      for (int f = 0; f < 2 * nd; f++)
        if (fc.plateOn[f] && ijk[f / 2] == ((f % 2) ? grid.n[f / 2] : 0)) vs[f / 2] = 0.0;
    }
    for (int c = 0; c < nd; c++) gvs[nd * i + c] = vs[c];
  }
  for (int f = 0; f < 6; f++) gvsU[f] = (minvU[f] != 0.0) ? gvsU[f] / mU[f] : 0.0;
  // velocity gradient, F, the kernel, the constitutive update
  tF.resize(P.F.size()); tR.resize(P.R.size()); tV.resize(N);
  teps.resize(P.eps.size()); tsh.resize(P.sigHat.size()); tsig.resize(P.sig.size());
  for (int p = 0; p < N; p++) {
    const double *Fn = &P.F[nn * p];
    double L[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, dF[9], F1[9];
    for (long k = sS[p]; k < sS[p + 1]; k++) {
      const int i = eIdx[k];
      for (int c = 0; c < nd; c++) {
        const double vv = (i >= 0) ? gvs[nd * i + c] : ((c == (-i - 1) / 2) ? gvsU[-i - 1] : 0.0);
        if (vv == 0.0) continue;
        const double *g = &sdSc[(nc * k + c) * nd];
        for (int j = 0; j < nd; j++) L[nd * c + j] += vv * g[j];
      }
    }
    for (int a = 0; a < nd; a++)
      for (int c = 0; c < nd; c++) dF[nd * a + c] = (a == c ? 1.0 : 0.0) + dt * L[nd * a + c];
    for (int a = 0; a < nd; a++)
      for (int c = 0; c < nd; c++) {
        double s = 0.0;
        for (int e = 0; e < nd; e++) s += dF[nd * a + e] * Fn[nd * e + c];
        F1[nd * a + c] = s;
      }
    double deps[6], R1[9], J1 = 0.0;
    int st;
    if (par.smallKin) {                 // ruling 9064acca3 (2): Delta eps = sym(dt L), V = V0
      st = smallIncrement(nd, dF, deps, R1) ? 0 : -1;
      J1 = 1.0;
    } else {
      st = (nd == 3) ? ladruno_hypo::hypoIncrementFromF(Fn, F1, deps, R1, &J1)
                     : ladruno_hypo::hypoIncrementFromFPlaneStrain(Fn, F1, deps, R1, &J1);
    }
    if (st != 0) {
      err = "the kinematic update of particle " + std::to_string((long long)P.id[p]) +
            " refused: det(F_mid) <= 0 or det(F) <= 0 (inverted or collapsed)";
      revertAll(mats, p);
      return -1;
    }
    double e1[6], sm[6];
    for (int k = 0; k < nvv; k++) e1[k] = P.eps[nvv * p + k] + deps[k];
    if (mats.setTrialStrain(p, e1, sm) < 0) {
      err = "the material of particle " + std::to_string((long long)P.id[p]) +
            " failed to integrate the trial strain";
      revertAll(mats, p + 1);
      return -1;
    }
    for (int k = 0; k < nvv; k++) {
      teps[nvv * p + k] = e1[k];
      tsh[nvv * p + k] = P.sig0[nvv * p + k] + sm[k];
    }
    if (par.smallKin) for (int k = 0; k < nvv; k++) tsig[nvv * p + k] = tsh[nvv * p + k];
    else if (nd == 3) ladruno_hypo::pushStress6(R1, &tsh[nvv * p], &tsig[nvv * p]);
    else ladruno_hypo::pushStressPlaneStrain(R1, &tsh[nvv * p], &tsig[nvv * p]);
    for (int k = 0; k < nn; k++) { tF[nn * p + k] = F1[k]; tR[nn * p + k] = R1[k]; }
    tV[p] = par.smallKin ? P.V0[p] : J1 * P.V0[p];
  }

  // ---- commit the particle materials; a refusal aborts (W17) -------------------------------
  mats.clearRefusals();
  for (int p = 0; p < N; p++) mats.commitState(p);
  const int nRef = mats.pendingRefusals();
  if (nRef > 0) {
    err = std::to_string(nRef) + " particle material point(s) REFUSED this commit; sibling "
          "points have already committed, so the run is aborted (ADR 232 W17)";
    return -2;
  }
  P.x.swap(tx); P.v.swap(tv); P.F.swap(tF); P.R.swap(tR); P.V.swap(tV);
  P.eps.swap(teps); P.sigHat.swap(tsh); P.sig.swap(tsig);
  P.pw.swap(tpw);
  // the grid pressure of this step, for `mpm info gridPressure` (T0 iii)
  {
    std::vector<std::pair<long, int>> ord;
    for (size_t a = 0; a < pnode.size(); a++) ord.push_back(std::make_pair(grid.publicId(nodes[pnode[a]]), (int)a));
    std::sort(ord.begin(), ord.end());
    u.gNode.clear(); u.gP.clear(); u.gPhat.clear();
    for (const auto &o : ord) {
      u.gNode.push_back(o.first);
      u.gP.push_back(pa[o.second]);
      u.gPhat.push_back(phat[o.second]);
    }
  }
  return 0;
}

} // namespace ladruno_mpm

#endif
