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

// Ladruno: WP-253 -- mode F's threaded step (Ladruno_implementation/253_mpm_performance.md §5).
// Included by LadrunoMPMCore.h only; OpenSees-free like the rest of the core (the host hands in
// the parallel-for as a function pointer, Engine::parFor).
//
// THE RULE: every result is bit-identical to the serial Engine::advance at every thread count,
// by construction, not by tolerance.
//   * Particle-local loops (supports, G2P, dF / the kernel / the push) are independent per
//     particle and write only that particle's slots.
//   * Grid accumulations (P2G mass, momentum, external and internal force with the anti-locking
//     corner term, the second projection) are a node-centric GATHER over each node's entries in
//     ascending (particle, entry) order -- the serial scatter's exact summation order -- through
//     a node -> entry transpose built deterministically (per-block histograms, blocks in particle
//     order). A node appears at most once in a particle's tensor-product support, and the
//     particle's cell corners lie in its support (J3), so the corner term of particle p is applied
//     right after p's support entry at that node, as the serial loop interleaves it. A corner
//     outside the support (never seen; J3 says impossible) sends the step to the serial path.
//   * Grid-node loops are independent per node. The serial-order bookkeeping that DOES depend on
//     order -- the slaved ghosts' fold order -- is kept in the serial first-touch order (a node's
//     first touch is its lowest entry). The touched list itself is in discovery order: only
//     per-node loops and the reset walk it.
//   * Reductions are order-free: integer counts, max, the lowest particle id / index.
//   * Every per-entry expression is the serial one, written in the same order of operations.
//   * Materials: threaded only with matThreaded (the host's Update && Form allowlist check);
//     otherwise parallel kinematics, then the material calls serially in particle order (the
//     exact serial call sequence, stopping at the first failing particle), then a parallel push.
//   * Failures name the LOWEST failing particle with the serial message; the materials touched
//     are reverted (the serial revert set, plus any later particle whose trial ran in parallel).

#ifndef LadrunoMPMThreads_h
#define LadrunoMPMThreads_h

#include <climits>
#include <utility>

namespace ladruno_mpm {

// particleSupport's entry count, without the tensor product
inline int supportCount(const Grid &g, Basis b, const double *x, const double *lp)
{
  AxisSupport ax[3];
  int c = 1;
  for (int d = 0; d < g.ndm; d++) {
    axisSupport(g, b, d, x[d], lp[d], ax[d]);
    c *= ax[d].count;
  }
  return c;
}

// particleSupport written at an offset: the same entries, values and order
inline void particleSupportAt(const Grid &g, Basis b, const double *x, const double *lp,
                              long *node, double *S, double *dS)
{
  const int nd = g.ndm;
  AxisSupport ax[3];
  for (int d = 0; d < nd; d++) axisSupport(g, b, d, x[d], lp[d], ax[d]);
  const int nz = (nd == 3) ? ax[2].count : 1;
  long e = 0;
  for (int kz = 0; kz < nz; kz++)
    for (int ky = 0; ky < ax[1].count; ky++)
      for (int kx = 0; kx < ax[0].count; kx++) {
        int ijk[3] = {ax[0].idx[kx], ax[1].idx[ky], (nd == 3) ? ax[2].idx[kz] : 0};
        const double sx = ax[0].S[kx], sy = ax[1].S[ky], sz = (nd == 3) ? ax[2].S[kz] : 1.0;
        const double gx = ax[0].dS[kx], gy = ax[1].dS[ky], gz = (nd == 3) ? ax[2].dS[kz] : 0.0;
        node[e] = g.extId(ijk);
        S[e] = sx * sy * sz;
        dS[nd * e]     = gx * sy * sz;
        dS[nd * e + 1] = sx * gy * sz;
        if (nd == 3) dS[nd * e + 2] = sx * sy * gz;
        e++;
      }
}

// The split path's quiet span (Engine::quietBegin / quietEnd, the host's): begin() BEFORE the
// fork of the region that precedes a serial phase, so that region's workers sleep at its join
// instead of spinning through the serial phase; end() before the next region's fork. Timing
// only: nothing it does reaches a result. The destructor ends an open span on every return.
struct QuietSpan {
  const Engine &e;
  int tok;
  bool on;
  explicit QuietSpan(const Engine &eng) : e(eng), tok(-1), on(false) {}
  ~QuietSpan() { end(); }
  void begin()
  {
    if (on || e.quietBegin == nullptr) return;
    tok = e.quietBegin();
    on = true;
  }
  void end()
  {
    if (!on) return;
    if (e.quietEnd != nullptr) e.quietEnd(tok);
    on = false;
  }
  QuietSpan(const QuietSpan &) = delete;
  QuietSpan &operator=(const QuietSpan &) = delete;
};

// f(lo, hi, b) over nb contiguous blocks of [0, n), block b = [n b / nb, n (b+1) / nb)
template <class F> void Engine::forBlocks(long n, int nb, const F &f) const
{
  if (n <= 0 || nb <= 0) return;
  struct Ctx { const F *f; long n; int nb; } c = {&f, n, nb};
  auto body = [](int b, void *v) -> int {
    const Ctx &cx = *static_cast<const Ctx *>(v);
    const long lo = (long)((long long)cx.n * b / cx.nb);
    const long hi = (long)((long long)cx.n * (b + 1) / cx.nb);
    (*cx.f)(lo, hi, b);
    return 0;
  };
  parFor(nb, nThreads, 1, body, &c, nullptr);
}

inline int Engine::blocksFor(long n) const
{
  return (int)std::min<long>(n, 8L * (long)nThreads);
}

inline int Engine::firstInadmissibleThreaded() const
{
  const int N = P.size(), nd = grid.ndm, nb = blocksFor(N);
  if (N == 0) return -1;
  std::vector<int> w((size_t)nb, -1);
  forBlocks(N, nb, [&](long lo, long hi, int b) {
    int worst = -1;
    for (int p = (int)lo; p < (int)hi; p++)
      for (int d = 0; d < nd; d++)
        if (!admissibleAxis(grid, par.basis, d, P.x[nd * p + d], P.lp[nd * p + d])) {
          if (worst < 0 || P.id[p] < P.id[worst]) worst = p;
          break;
        }
    w[b] = worst;
  });
  int worst = -1;
  for (int b = 0; b < nb; b++)
    if (w[b] >= 0 && (worst < 0 || P.id[w[b]] < P.id[worst])) worst = w[b];
  return worst;
}

inline double Engine::criticalTimeStepThreaded(MaterialPoints &mats) const
{
  const int N = P.size(), nd = grid.ndm, nb = blocksFor(N);
  std::vector<double> bm((size_t)std::max(nb, 1), 0.0);
  forBlocks(N, nb, [&](long lo, long hi, int b) {
    double wmax = 0.0;
    for (int p = (int)lo; p < (int)hi; p++) {
      const double rho = P.m[p] / P.V[p];
      const double M = mats.normalModulus(p);
      const double c = (M > 0.0) ? std::sqrt(M / rho) : 0.0;
      double v2 = 0.0;
      for (int d = 0; d < nd; d++) v2 += P.v[nd * p + d] * P.v[nd * p + d];
      const double w = c + std::sqrt(v2);
      if (w > wmax) wmax = w;
    }
    bm[b] = wmax;
  });
  double wmax = 0.0;
  for (int b = 0; b < nb; b++)
    if (bm[b] > wmax) wmax = bm[b];
  if (!(wmax > 0.0)) return std::numeric_limits<double>::infinity();
  return grid.h / wmax;
}

inline double Engine::massRefThreaded() const
{
  const int N = P.size(), nb = blocksFor(N);
  std::vector<double> bm((size_t)std::max(nb, 1), 0.0);
  forBlocks(N, nb, [&](long lo, long hi, int b) {
    double r = 0.0;
    for (long p = lo; p < hi; p++) r = std::max(r, P.m[p]);
    bm[b] = r;
  });
  double r = 0.0;
  for (int b = 0; b < nb; b++) r = std::max(r, bm[b]);
  return r;
}

// buildSupports, threaded: count, prefix (per block), fill at the offsets. Identical arrays.
inline void Engine::buildSupportsThreaded()
{
  const int N = P.size(), nd = grid.ndm, nb = blocksFor(N);
  supCnt.resize((size_t)N);
  tBlkSum.assign((size_t)nb + 1, 0);
  forBlocks(N, nb, [&](long lo, long hi, int b) {
    long s = 0;
    for (long p = lo; p < hi; p++) {
      supCnt[p] = supportCount(grid, par.basis, &P.x[nd * p], &P.lp[nd * p]);
      s += supCnt[p];
    }
    tBlkSum[b + 1] = s;
  });
  for (int b = 0; b < nb; b++) tBlkSum[b + 1] += tBlkSum[b];
  const long total = tBlkSum[nb];
  supStart.resize((size_t)N + 1);
  supStart[0] = 0;
  supNode.resize((size_t)total);
  supS.resize((size_t)total);
  supdS.resize((size_t)nd * total);
  pOf.resize((size_t)total);
  forBlocks(N, nb, [&](long lo, long hi, int b) {
    long off = tBlkSum[b];
    for (long p = lo; p < hi; p++) {
      particleSupportAt(grid, par.basis, &P.x[nd * p], &P.lp[nd * p], &supNode[off], &supS[off],
                        &supdS[nd * off]);
      for (int e = 0; e < supCnt[p]; e++) pOf[off + e] = (int)p;
      off += supCnt[p];
      supStart[p + 1] = off;
    }
  });
  // the grid arrays, exactly as buildSupports sizes them
  const long nN = grid.numNodes();
  if ((long)gm.size() != nN) {
    gm.assign(nN, 0.0);
    gp.assign(nd * nN, 0.0); gfi.assign(nd * nN, 0.0); gfe.assign(nd * nN, 0.0);
    gv1.assign(nd * nN, 0.0); ga.assign(nd * nN, 0.0); gvs.assign(nd * nN, 0.0);
    gmark.assign(nN, 0); gact.assign(nN, 0); gslave.assign(nN, 0);
    touched.clear();
  }
  if (tMark.n != nN) {
    tMark.a.reset(new std::atomic<int>[(size_t)nN]);
    for (long i = 0; i < nN; i++) tMark.a[i].store(0, std::memory_order_relaxed);
    tMark.n = nN;
  }
  if ((long)tSlot.size() != nN) tSlot.assign((size_t)nN, 0);
}

// The touched nodes (discovery order) and the node -> entry transpose in ascending entry order.
// One block per thread. A block's histogram row is WINDOWED (WP-253 review): it spans only the
// slots of the owner blocks (the discovering blocks) of the nodes it touches, not all nTouched
// slots. Particles that are spatially coherent give windows of about 3 nTouched / nb, so memory
// and the per-thread clear/sum work scale with 1/threads; incoherent particles widen every window
// to all owners, which is the dense nb x nTouched histogram this replaced and is never worse.
// Slots, tPtr and tEntry are the same as the dense version's: only the histogram storage changed.
inline void Engine::transposeThreaded()
{
  const int N = P.size();
  const int nbT = std::max(1, std::min(N, nThreads));
  if ((int)tDisc.size() < nbT) tDisc.resize((size_t)nbT);
  tWin.resize((size_t)2 * nbT);
  // (a) discovery: each node is claimed by exactly one block, which stores (block + 1) in its
  //     mark; every block records the range of owner blocks among the nodes it touches
  forBlocks(N, nbT, [&](long lo, long hi, int b) {
    std::vector<long> &d = tDisc[b];
    d.clear();
    int oLo = b, oHi = b;
    for (long k = supStart[lo]; k < supStart[hi]; k++) {
      const long i = supNode[k];
      // a compare-exchange, never an exchange: a losing claim must not overwrite the winner's
      // id (a later reader would take the wrong owner, and its window would miss the node)
      int o = tMark.a[i].load(std::memory_order_relaxed);
      if (o == 0 && tMark.a[i].compare_exchange_strong(o, b + 1, std::memory_order_relaxed)) {
        d.push_back(i);
        gmark[i] = 1;
        o = b + 1;
      }
      o -= 1;
      if (o < oLo) oLo = o;
      if (o > oHi) oHi = o;
    }
    tWin[2 * b] = oLo;
    tWin[2 * b + 1] = oHi;
  });
  // GUARD (WP-253): `touched` is in DISCOVERY order, which depends on thread timing (which block
  // claims a node first) and is NOT the serial first-touch order. Every loop that walks it must
  // therefore be per-node (no cross-node accumulation, no order-dependent bookkeeping). The one
  // order-dependent consumer, K4's slaved-ghost list, re-sorts by each node's first entry
  // (tEntry[tPtr[t]], = the serial first touch) before folding. A new walker that needs order
  // must do the same.
  touched.clear();
  for (int b = 0; b < nbT; b++) touched.insert(touched.end(), tDisc[b].begin(), tDisc[b].end());
  const long nT = (long)touched.size();
  nTouchedP2G = (int)nT;
  for (long t = 0; t < nT; t++) tSlot[touched[t]] = (int)t;
  // owner o's slots are [tOwnOff[o], tOwnOff[o+1]); block b's row covers its owner window and
  // row(b)[t - tOwnOff[tWin[2b]]] is its count (later its write cursor) for slot t
  tOwnOff.resize((size_t)nbT + 1);
  tOwnOff[0] = 0;
  for (int o = 0; o < nbT; o++) tOwnOff[o + 1] = tOwnOff[o] + (long)tDisc[o].size();
  tHOff.resize((size_t)nbT + 1);
  tHOff[0] = 0;
  for (int b = 0; b < nbT; b++)
    tHOff[b + 1] = tHOff[b] + (tOwnOff[tWin[2 * b + 1] + 1] - tOwnOff[tWin[2 * b]]);
  tCovPtr.assign((size_t)nbT + 1, 0);
  tCov.clear();
  for (int o = 0; o < nbT; o++) {
    for (int b = 0; b < nbT; b++)
      if (tWin[2 * b] <= o && o <= tWin[2 * b + 1]) tCov.push_back(b);
    tCovPtr[o + 1] = (int)tCov.size();
  }
  auto cell = [&](int b, long t) -> int & {
    return tHist[(size_t)(tHOff[b] + t - tOwnOff[tWin[2 * b]])];
  };
  // (b) per-block histograms over the block's window
  tHist.resize((size_t)tHOff[nbT]);
  forBlocks(N, nbT, [&](long lo, long hi, int b) {
    std::fill(tHist.begin() + tHOff[b], tHist.begin() + tHOff[b + 1], 0);
    for (long k = supStart[lo]; k < supStart[hi]; k++) cell(b, tSlot[supNode[k]])++;
  });
  // (c) node totals, their prefix, then each block's write offset per node (blocks ascending:
  //     the blocks whose window holds the node's owner, every other block's count being zero)
  tPtr.resize((size_t)nT + 1);
  tPtr[0] = 0;
  const int nbN = blocksFor(nT);
  auto ownerOf = [&](long t) -> int {   // the owner block of slot t (t < nT)
    return (int)(std::upper_bound(tOwnOff.begin(), tOwnOff.end(), t) - tOwnOff.begin()) - 1;
  };
  forBlocks(nT, nbN, [&](long lo, long hi, int) {
    if (lo >= hi) return;
    int o = ownerOf(lo);
    for (long t = lo; t < hi; t++) {
      while (t >= tOwnOff[o + 1]) o++;
      int s = 0;
      for (int c = tCovPtr[o]; c < tCovPtr[o + 1]; c++) s += cell(tCov[c], t);
      tPtr[t + 1] = s;
    }
  });
  for (long t = 0; t < nT; t++) tPtr[t + 1] += tPtr[t];
  forBlocks(nT, nbN, [&](long lo, long hi, int) {
    if (lo >= hi) return;
    int o = ownerOf(lo);
    for (long t = lo; t < hi; t++) {
      while (t >= tOwnOff[o + 1]) o++;
      int run = tPtr[t];
      for (int c = tCovPtr[o]; c < tCovPtr[o + 1]; c++) {
        int &h = cell(tCov[c], t);
        const int n = h;
        h = run;
        run += n;
      }
    }
  });
  // (d) fill: a block writes its entries in ascending k into its own runs
  tEntry.resize((size_t)supStart[N]);
  forBlocks(N, nbT, [&](long lo, long hi, int b) {
    for (long k = supStart[lo]; k < supStart[hi]; k++) tEntry[cell(b, tSlot[supNode[k]])++] = (int)k;
  });
}

inline void Engine::clearGridThreaded()
{
  const int nd = grid.ndm;
  const long nT = (long)touched.size();
  forBlocks(nT, blocksFor(nT), [&](long lo, long hi, int) {
    for (long t = lo; t < hi; t++) {
      const long i = touched[t];
      gm[i] = 0.0;
      gmark[i] = 0;
      gact[i] = 0;
      gslave[i] = 0;
      for (int a = 0; a < nd; a++)
        gp[nd * i + a] = gfi[nd * i + a] = gfe[nd * i + a] =
        gv1[nd * i + a] = ga[nd * i + a] = gvs[nd * i + a] = 0.0;
      tMark.a[i].store(0, std::memory_order_relaxed);
    }
  });
  touched.clear();
  slave.clear(); slaveImg.clear(); slaveSg.clear();
  nTouchedP2G = 0;
}

// The serial advance for this one step (nothing of the threaded step has touched any state the
// serial path reads: the caller has reset the grid). lastXferThreads reports 1.
inline int Engine::advanceSerialFallback(double dt, double lambda, MaterialPoints &mats,
                                         std::string &err)
{
  const int nt = nThreads;
  nThreads = 1;
  const int r = advance(dt, lambda, mats, err);
  nThreads = nt;
  return r;
}

inline int Engine::advanceThreaded(double dt, double lambda, MaterialPoints &mats, std::string &err)
{
  const int N = P.size(), nd = grid.ndm, nvv = nv();
  lastMatFail = -1;
  if (N == 0) return 0;
  if (!bodyForce.empty() && (long)bodyForce.size() != (long)nd * N) {
    err = "the particle body force holds " + std::to_string((long long)bodyForce.size()) +
          " values for " + std::to_string(N) + " particles (ndm per particle): set it again;"
          " step refused, nothing changed";
    return -1;
  }
  lastXferThreads = nThreads;
  lastMatThreads = matThreaded ? nThreads : 1;
  const int nbP = blocksFor(N);
  PhaseRun ph(phaseHook);             // WP-253: profiler phases (timing only), calling thread
  QuietSpan quiet(*this);             // the split path's sleeping workers (timing only)
  ph.go(PH_SUPPORT);
  buildSupportsThreaded();
  if (supStart[N] >= (long)INT_MAX) {  // the transpose indexes entries with int
    ph.stop();
    return advanceSerialFallback(dt, lambda, mats, err);
  }

  // ---- step 3: P2G, as a gather ----------------------------------------------
  ph.go(PH_TRANSPOSE);
  transposeThreaded();
  ph.go(PH_P2G);
  const long nT = nTouchedP2G;
  const double mref = massRefThreaded();          // == massRef(): a max is order-free
  double b[3];
  for (int d = 0; d < 3; d++) b[d] = lambda * par.gravity[d];
  const bool al = par.antiLocking != 0;
  auto massNode = [&](long t) {
    const long i = touched[t];
    for (int e = tPtr[t]; e < tPtr[t + 1]; e++) {
      const long k = tEntry[e];
      gm[i] += supS[k] * P.m[pOf[k]];
    }
  };
  auto forceNode = [&](long t) {
    const long i = touched[t];
    for (int e = tPtr[t]; e < tPtr[t + 1]; e++) {
      const long k = tEntry[e];
      const int p = pOf[k];
      const double mp = P.m[p], Vp = P.V[p];
      const double *bfp = bodyForce.empty() ? nullptr : &bodyForce[nd * p];
      double sg[3][3];
      voigtToTensor(&P.sig[nvv * p], sg);
      const bool avg = par.antiLocking && talFlag[p] == 1;
      const double pbar = avg ? pbarP[p] : 0.0;
      const double S = supS[k];
      const double *g = &supdS[nd * k];
      for (int a = 0; a < nd; a++) {
        gp[nd * i + a]   += S * mp * P.v[nd * p + a];
        gfe[nd * i + a]  += S * mp * (bfp ? b[a] + bfp[a] : b[a]);
        if (LADRUNO_MUTATE_MPM != MUT_ZERO) {
          double s = 0.0;
          for (int c = 0; c < nd; c++) s += sg[a][c] * g[c];
          if (avg) s -= pbar * g[a];
          gfi[nd * i + a] -= Vp * s;
        }
      }
      // the corner term of p at this node, right after p's support entry here (serial order)
      if (avg && LADRUNO_MUTATE_MPM != MUT_ZERO && cornerOfEntry[k] >= 0) {
        const int c = cornerOfEntry[k];
        const double ac = (nd == 3) ? 0.25 / grid.h : 0.5 / grid.h;   // cellCorners' gN
        for (int a = 0; a < nd; a++) {
          const double gNa = ((c >> a) & 1) ? ac : -ac;
          gfi[nd * i + a] -= Vp * pbar * gNa;
        }
      }
    }
  };
  const int nbN = blocksFor(nT);
  if (!al) {
    forBlocks(nT, nbN, [&](long lo, long hi, int) {
      for (long t = lo; t < hi; t++) { massNode(t); forceNode(t); }
    });
  } else {
    forBlocks(nT, nbN, [&](long lo, long hi, int) {
      for (long t = lo; t < hi; t++) massNode(t);
    });
    // WP-232j (J5): each particle's path from the P2G mass, its pbar, and where its corners
    // sit among its support entries
    const double mcut0 = par.masscut * mref;
    talFlag.resize((size_t)N);
    pbarP.resize((size_t)N);
    cornerOfEntry.resize((size_t)supStart[N]);
    blkFirst.assign((size_t)nbP, 0);
    forBlocks(N, nbP, [&](long lo, long hi, int bl) {
      for (long p = lo; p < hi; p++) {
        long ce[8];
        double gN[8][3];
        const int nc = cellCorners(grid, &P.x[nd * p], ce, gN);
        unsigned char f = 1;
        for (int k = 0; k < nc; k++)
          if (!(gm[ce[k]] > 0.0 && imageCutoffMass(grid, par.basis, fc, ce[k], gm[ce[k]]) > mcut0)) {
            f = 2;
            break;
          }
        talFlag[p] = f;
        for (long k = supStart[p]; k < supStart[p + 1]; k++) cornerOfEntry[k] = -1;
        if (f != 1) continue;
        double sg[3][3];
        voigtToTensor(&P.sig[nvv * p], sg);
        double pbar = 0.0;
        for (int a = 0; a < nd; a++) pbar += sg[a][a];
        pbar /= nd;
        pbarP[p] = pbar;
        for (int c = 0; c < nc; c++) {
          long k = supStart[p];
          while (k < supStart[p + 1] && supNode[k] != ce[c]) k++;
          if (k == supStart[p + 1]) { blkFirst[bl] = 1; break; }   // J3 broken: serial step
          cornerOfEntry[k] = (signed char)c;
        }
      }
    });
    for (int bl = 0; bl < nbP; bl++)
      if (blkFirst[bl]) {
        ph.stop();
        clearGridThreaded();
        return advanceSerialFallback(dt, lambda, mats, err);
      }
    forBlocks(nT, nbN, [&](long lo, long hi, int) {
      for (long t = lo; t < hi; t++) forceNode(t);
    });
  }

  // K4 (B6): the slaved ghosts in the serial first-touch order (= ascending first entry), then
  // the folds, serially, exactly as collectImages and the fold loop
  if (isBSpline(par.basis)) {
    tImg.resize((size_t)nT);
    tSg.resize((size_t)3 * nT);
    forBlocks(nT, nbN, [&](long lo, long hi, int) {
      for (long t = lo; t < hi; t++) {
        long img = -1;
        tImg[t] = imageOfNode(grid, par.basis, fc, touched[t], img, &tSg[3 * t]) ? img : -1;
      }
    });
    std::vector<std::pair<int, long>> order;
    for (long t = 0; t < nT; t++)
      if (tImg[t] >= 0) order.push_back(std::make_pair(tEntry[tPtr[t]], t));
    std::sort(order.begin(), order.end());
    for (size_t q = 0; q < order.size(); q++) {
      const long t = order[q].second, i = touched[t], img = tImg[t];
      gslave[i] = 1;
      slave.push_back(i);
      slaveImg.push_back(img);
      for (int d = 0; d < 3; d++) slaveSg.push_back(tSg[3 * t + d]);
      if (!gmark[img]) { gmark[img] = 1; touched.push_back(img); }
    }
    for (size_t s = 0; s < slave.size(); s++) {
      gm[slaveImg[s]] += gm[slave[s]];
      foldVec(gp, s);
      foldVec(gfi, s);
      foldVec(gfe, s);
    }
    // zeroFaceOdd(gp / gfi / gfe), per node
    const long nAll = (long)touched.size();
    forBlocks(nAll, blocksFor(nAll), [&](long lo, long hi, int) {
      for (long t = lo; t < hi; t++) {
        const long i = touched[t];
        int ijk[3];
        grid.extIJK(i, ijk);
        for (int d = 0; d < nd; d++) {
          unsigned m = 0u;
          if (ijk[d] == 0) m |= fc.fixMask[2 * d];
          if (ijk[d] == grid.n[d]) m |= fc.fixMask[2 * d + 1];
          for (int a = 0; a < nd; a++)
            if (m & (1u << a)) gp[nd * i + a] = gfi[nd * i + a] = gfe[nd * i + a] = 0.0;
        }
      }
    });
  }
  const long nAll = (long)touched.size();
  const int nbA = blocksFor(nAll);

  // ---- step 4: damping, grid momentum update, grid conditions (per node) -------
  ph.go(PH_GRID);
  const double mcut = par.masscut * mref;
  forBlocks(nAll, nbA, [&](long lo, long hi, int) {
    for (long t = lo; t < hi; t++) {
      const long i = touched[t];
      if (gslave[i]) continue;          // K4: set from the image below
      const bool act = gm[i] > 0.0 && imageCutoffMass(grid, par.basis, fc, i, gm[i]) > mcut;
      gact[i] = act ? 1 : 0;
      if (!act) {
        for (int a = 0; a < nd; a++) gv1[nd * i + a] = ga[nd * i + a] = 0.0;
        continue;
      }
      double vn[3], v1[3];
      for (int a = 0; a < nd; a++) {
        vn[a] = gp[nd * i + a] / gm[i];
        double f = gfi[nd * i + a] + gfe[nd * i + a];
        const double sg = (vn[a] > 0.0) ? 1.0 : ((vn[a] < 0.0) ? -1.0 : 0.0);
        f -= par.damping * std::fabs(f) * sg;
        v1[a] = vn[a] + dt * f / gm[i];
      }
      applyConditions(i, v1, true);
      for (int a = 0; a < nd; a++) {
        gv1[nd * i + a] = v1[a];
        ga[nd * i + a] = (v1[a] - vn[a]) / dt;
      }
    }
  });
  for (size_t s = 0; s < slave.size(); s++) {
    gact[slave[s]] = gact[slaveImg[s]];
    reflectVec(gv1, s);
    reflectVec(ga, s);
  }

  // ---- step 5: MUSL particle update (trial), with A7's check (per particle) -----
  ph.go(PH_G2P);
  tx.resize(P.x.size());
  tv.resize(P.v.size());
  blkFirst.assign((size_t)nbP, -1);
  forBlocks(N, nbP, [&](long lo, long hi, int bl) {
    int worst = -1;
    for (int p = (int)lo; p < (int)hi; p++) {
      double sa[3] = {0, 0, 0}, sv[3] = {0, 0, 0};
      for (long k = supStart[p]; k < supStart[p + 1]; k++) {
        const long i = supNode[k];
        for (int a = 0; a < nd; a++) {
          sa[a] += supS[k] * ga[nd * i + a];
          sv[a] += supS[k] * gv1[nd * i + a];
        }
      }
      for (int a = 0; a < nd; a++) {
        const double vp = P.v[nd * p + a];
        tv[nd * p + a] = par.flip * (vp + dt * sa[a]) + (1.0 - par.flip) * sv[a];
        double dx = dt * sv[a];
        if (LADRUNO_MUTATE_MPM == MUT_SCALE) dx *= MUT_SCALE_FACTOR;
        tx[nd * p + a] = P.x[nd * p + a] + dx;
      }
      for (int a = 0; a < nd; a++)
        if (!admissibleAxis(grid, par.basis, a, tx[nd * p + a], P.lp[nd * p + a])) {
          if (worst < 0 || P.id[p] < P.id[worst]) worst = p;
          break;
        }
    }
    blkFirst[bl] = worst;
  });
  {
    int worst = -1;
    for (int bl = 0; bl < nbP; bl++)
      if (blkFirst[bl] >= 0 && (worst < 0 || P.id[blkFirst[bl]] < P.id[worst])) worst = blkFirst[bl];
    if (worst >= 0) {
      err = "the step would carry particle " + std::to_string((long long)P.id[worst]) +
            " out of the grid (interface A7); step refused, nothing changed";
      ph.stop();
      clearGridThreaded();
      return -1;
    }
  }

  // second projection of the updated particle momentum (gather), then its grid conditions
  ph.go(PH_PROJECT2);
  forBlocks(nAll, nbA, [&](long lo, long hi, int) {
    for (long t = lo; t < hi; t++) {
      const long i = touched[t];
      for (int a = 0; a < nd; a++) gvs[nd * i + a] = 0.0;
      if (t >= nT) continue;            // an image added by K4: no entries of its own
      for (int e = tPtr[t]; e < tPtr[t + 1]; e++) {
        const long k = tEntry[e];
        const int p = pOf[k];
        for (int a = 0; a < nd; a++)
          gvs[nd * i + a] += supS[k] * P.m[p] * tv[nd * p + a];
      }
    }
  });
  for (size_t s = 0; s < slave.size(); s++) foldVec(gvs, s);   // K4 (B6)
  forBlocks(nAll, nbA, [&](long lo, long hi, int) {
    for (long t = lo; t < hi; t++) {
      const long i = touched[t];
      if (gslave[i]) continue;
      double vs[3] = {0, 0, 0};
      if (gact[i]) {
        for (int a = 0; a < nd; a++) vs[a] = gvs[nd * i + a] / gm[i];
        applyConditions(i, vs, false);
      }
      for (int a = 0; a < nd; a++) gvs[nd * i + a] = vs[a];
    }
  });
  for (size_t s = 0; s < slave.size(); s++) reflectVec(gvs, s);

  // ---- velocity gradient, F, the kernel, the constitutive update --------------
  ph.go(PH_UPDATE);
  const int nn = nd * nd;
  tF.resize(P.F.size()); tR.resize(P.R.size()); tV.resize(N);
  teps.resize(P.eps.size()); tsh.resize(P.sigHat.size()); tsig.resize(P.sig.size());
  upStatus.resize((size_t)N);
  // the kinematics of particle p: 0, 1 (anti-locking refused) or 2 (kernel refused); on 0 it
  // leaves e1 in teps and F1, R1, J1 V0 in tF, tR, tV -- the serial loop's values
  auto kinematics = [&](int p) -> unsigned char {
    const double *Fn = &P.F[nn * p];
    double dF[9], dFu[9], F1[9];
    if (stepIncrement(p, gvs, dt, par.antiLocking ? talFlag[p] : 0, dF, dFu) != 0) return 1;
    for (int a = 0; a < nd; a++)
      for (int c = 0; c < nd; c++) {
        double s = 0.0;
        for (int e = 0; e < nd; e++) s += dF[nd * a + e] * Fn[nd * e + c];
        F1[nd * a + c] = s;
      }
    double deps[6], R1[9], J1 = 0.0;
    int st;
    if (par.smallKin) {              // WP-232l (P6): small strain, R = I, V = V0
      st = smallIncrement(nd, dF, deps, R1) ? 0 : -1;
      J1 = 1.0;
    } else {
      st = (nd == 3) ? ladruno_hypo::hypoIncrementFromF(Fn, F1, deps, R1, &J1)
                     : ladruno_hypo::hypoIncrementFromFPlaneStrain(Fn, F1, deps, R1, &J1);
    }
    if (st != 0) return 2;
    for (int k = 0; k < nvv; k++) teps[nvv * p + k] = P.eps[nvv * p + k] + deps[k];
    for (int k = 0; k < nn; k++) { tF[nn * p + k] = F1[k]; tR[nn * p + k] = R1[k]; }
    tV[p] = J1 * P.V0[p];
    return 0;
  };
  auto push = [&](int p) {
    if (par.smallKin) for (int k = 0; k < nvv; k++) tsig[nvv * p + k] = tsh[nvv * p + k];
    else if (nd == 3) ladruno_hypo::pushStress6(&tR[nn * p], &tsh[nvv * p], &tsig[nvv * p]);
    else ladruno_hypo::pushStressPlaneStrain(&tR[nn * p], &tsh[nvv * p], &tsig[nvv * p]);
  };
  auto failText = [&](int p, int kind) {
    const std::string id = std::to_string((long long)P.id[p]);
    if (kind == 1)
      return "the anti-locking increment of particle " + id +
             " refused: det(dF) <= 0 or det(dF_cell) <= 0 (inverted or collapsed)";
    if (kind == 2)
      return "the kinematic update of particle " + id +
             " refused: det(F_mid) <= 0 or det(F) <= 0 (inverted or collapsed)";
    return "the material of particle " + id + " failed to integrate the trial strain";
  };
  if (matThreaded) {
    // every material on the allowlist: kinematics, material and push in one particle loop
    blkFirst.assign((size_t)nbP, -1);
    forBlocks(N, nbP, [&](long lo, long hi, int bl) {
      int first = -1;
      for (int p = (int)lo; p < (int)hi; p++) {
        unsigned char st = kinematics(p);
        if (st == 0) {
          double sm[6];
          if (mats.setTrialStrain(p, &teps[nvv * p], sm) < 0) st = 3;
          else {
            for (int k = 0; k < nvv; k++) tsh[nvv * p + k] = P.sig0[nvv * p + k] + sm[k];
            push(p);
          }
        }
        upStatus[p] = st;
        if (st != 0 && first < 0) first = p;
      }
      blkFirst[bl] = first;
    });
    mats.afterParallel();
    int pf = -1;
    for (int bl = 0; bl < nbP && pf < 0; bl++) pf = blkFirst[bl];
    if (pf >= 0) {
      const int kind = upStatus[pf];
      err = failText(pf, kind);
      if (kind == 3) lastMatFail = pf;
      // the serial revert set (q < pf, and pf itself after its material ran), plus every later
      // particle whose trial ran in parallel (status 0 or 3)
      for (int q = 0; q < N; q++)
        if (q < pf || (q == pf && kind == 3) ||
            (q > pf && (upStatus[q] == 0 || upStatus[q] == 3)))
          mats.revertToLastCommit(q);
      ph.stop();
      clearGridThreaded();
      return -1;
    }
  } else {
    // parallel kinematics, the material calls serially in particle order (the serial call
    // sequence, stopping at the first failing particle), then a parallel push. The kinematics
    // region forks inside a quiet span, so its workers sleep through the serial calls.
    quiet.begin();
    forBlocks(N, nbP, [&](long lo, long hi, int) {
      for (int p = (int)lo; p < (int)hi; p++) upStatus[p] = kinematics(p);
    });
    for (int p = 0; p < N; p++) {
      if (upStatus[p] != 0) {
        err = failText(p, upStatus[p]);
        revertAll(mats, p);
        ph.stop();
        clearGridThreaded();
        return -1;
      }
      double sm[6];
      if (mats.setTrialStrain(p, &teps[nvv * p], sm) < 0) {
        err = failText(p, 3);
        revertAll(mats, p + 1);
        ph.stop();
        clearGridThreaded();
        return -1;
      }
      for (int k = 0; k < nvv; k++) tsh[nvv * p + k] = P.sig0[nvv * p + k] + sm[k];
    }
    quiet.end();                      // the push and later regions keep the previous wait
    forBlocks(N, nbP, [&](long lo, long hi, int) {
      for (int p = (int)lo; p < (int)hi; p++) push(p);
    });
  }
  ph.stop();                          // the grid reset is the advance's own time
  if (!matThreaded) quiet.begin();    // the reset's workers sleep through the serial commit
  clearGridThreaded();

  // ---- step 6: commit the particle materials; a refusal aborts (W17) -----
  ph.go(PH_COMMIT);
  mats.clearRefusals();
  if (matThreaded) {
    forBlocks(N, nbP, [&](long lo, long hi, int) {
      for (int p = (int)lo; p < (int)hi; p++) mats.commitState(p);
    });
    mats.afterParallel();
  } else {
    for (int p = 0; p < N; p++) mats.commitState(p);
  }
  quiet.end();
  const int nRef = mats.pendingRefusals();
  ph.stop();
  if (nRef > 0) {
    err = std::to_string(nRef) + " particle material point(s) REFUSED this commit; sibling "
          "points have already committed, so the run is aborted (ADR 232 W17)";
    return -2;
  }

  // ---- step 7: G2P kinematics take their trial values ---------------------
  P.x.swap(tx); P.v.swap(tv); P.F.swap(tF); P.R.swap(tR); P.V.swap(tV);
  P.eps.swap(teps); P.sigHat.swap(tsh); P.sig.swap(tsig);
  if (par.antiLocking) countAntiLocking(talFlag);
  return 0;
}

} // namespace ladruno_mpm

#endif
