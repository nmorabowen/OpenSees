// Standalone driver for the shared Biot u-p (LadrunoUP) assembly kernel — NO
// OpenSees deps. Reimplements single-element block assembly using ONLY the pure
// headers (SRC/element/ladrunoUP/LadrunoUPKernel.h + LadrunoUPShapes.h). It:
//
//   (1) ingests an oracle-emitted cases file (argv[1], default
//       tests/ladruno_up_cases.txt), drives the matching shape provider over its
//       PINNED GP rule, folds dv = w·detJ·thick, accumulates Q/H/S/H̃/f_seep via
//       the kernel API into zeroed buffers, runs dofMap, and compares EVERYTHING
//       against the case's stored values to ≤1e-9 (worst delta per block per
//       case; nonzero exit on any failure);
//   (2) always runs self-tests independent of the cases file — provider
//       invariants (partition of unity, Σ∇N = 0, Σ w·detJ = volume, detJ > 0,
//       TH pressure partition) plus the scalar kernel helpers stabAlphaAuto /
//       oneOverQbar / elemSizeH against hand-computed contract-formula values.
//
// The companion pytest (test_ladruno_up_kernel_cpp.py) compiles + runs this with
// g++ and asserts exit 0 — catching any C++ transcription bug in the coupled
// u-p block math without building OpenSees (the ADR-71 P0 novel-math surface).
// Independence: this file NEVER reads the numpy oracle (ladruno_up_reference.py);
// two independent implementations agreeing IS the verification (WP0.D, MAIN).
//
// CONTRACT NOTE (flagged for MAIN): the plan addendum §5 folds dv = w·detJ·thick
// with the literal (signed) detJ the provider returns. For all canonical
// (properly-oriented, detJ > 0) cases this equals w·|detJ|·thick; a negatively-
// wound element would sign-flip every block identically on both sides, so the
// cross-check still holds. Self-test geometries are positively oriented.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>

#include <LadrunoUPKernel.h>
#include <LadrunoUPShapes.h>

using namespace ladruno_up;

static const double BLOCK_TOL = 1e-9;   // mixed abs/rel gate on stored blocks
static int g_fail = 0;

// ------------------------------------------------------------------ parsing -- //
struct Matrix {
  int rows = 0, cols = 0;
  std::vector<double> v;
  double at(int r, int c) const { return v[(size_t)r * cols + c]; }
};
struct Case {
  std::string name, shape, porder;
  int ndm = 0, nNodes = 0;
  std::vector<double> coords;               // nNodes*ndm
  double biotAlpha = 0, oneOverQbar = 0, rhoF = 0, stabAlpha = 0, thick = 1;
  std::vector<double> kbar, drive;          // ndm each
  int totalDof = 0;
  std::vector<int> uOff, pOff;              // nNodes each
  Matrix Q, H, S, HT;
  std::vector<double> FSEEP;
};

static Matrix readMat(std::istream& in) {
  // rows cols then rows*cols floats (token stream)
  Matrix m;
  in >> m.rows >> m.cols;
  m.v.resize((size_t)m.rows * m.cols);
  for (auto& x : m.v) in >> x;
  return m;
}

// Report the worst mismatch of a computed buffer vs a stored one.
static void compare(const char* tag, const Case& c, const std::vector<double>& got,
                    const std::vector<double>& want) {
  if (got.size() != want.size()) {
    printf("  [%s] %s SIZE MISMATCH got=%zu want=%zu\n", c.name.c_str(), tag,
           got.size(), want.size());
    g_fail++;
    return;
  }
  double worstAbs = 0.0, worstMetric = 0.0;
  for (size_t i = 0; i < got.size(); i++) {
    double a = std::fabs(got[i] - want[i]);
    double metric = a / (1.0 + std::fabs(want[i]));
    if (a > worstAbs) worstAbs = a;
    if (metric > worstMetric) worstMetric = metric;
  }
  const char* verdict = (worstMetric <= BLOCK_TOL) ? "ok" : "FAIL";
  if (worstMetric > BLOCK_TOL) g_fail++;
  printf("  [%s] %-6s worstAbs=%.3e worstMetric=%.3e  %s\n",
         c.name.c_str(), tag, worstAbs, worstMetric, verdict);
}

// ----------------------------------------------------------- provider driver -- //
template <class P>
static void driveCase(const Case& c) {
  if (c.ndm != P::ndm || c.nNodes != P::nNodes) {
    printf("  [%s] shape/dim mismatch\n", c.name.c_str());
    g_fail++;
    return;
  }
  const int ndm = P::ndm, nNodes = P::nNodes;
  const bool th = (c.porder == "th");

  int pCarrier[16];
  for (int n = 0; n < nNodes; n++) pCarrier[n] = th ? P::pCarrierTH[n] : 1;
  int nNp = 0;
  for (int n = 0; n < nNodes; n++) nNp += pCarrier[n];

  const int nu = nNodes * ndm;
  std::vector<double> Q((size_t)nu * nNp, 0.0), H((size_t)nNp * nNp, 0.0),
      S((size_t)nNp * nNp, 0.0), HT((size_t)nNp * nNp, 0.0), FS(nNp, 0.0);

  for (int gp = 0; gp < P::nGP; gp++) {
    double Nu[16], dNu[48], Np[16], dNp[48], detJ = 0.0;
    P::eval(gp, c.coords.data(), Nu, dNu, &detJ);
    P::evalP(gp, c.coords.data(), th, Np, dNp);
    double dv = P::gpWeight[gp] * detJ * c.thick;   // signed detJ, per addendum §5
    addQ(dNu, Np, nNodes, nNp, ndm, c.biotAlpha, dv, Q.data(), nNp);
    addH(dNp, nNp, ndm, c.kbar.data(), dv, H.data(), nNp);
    addS(Np, nNp, c.oneOverQbar, dv, S.data(), nNp);
    addHtilde(dNp, nNp, ndm, c.stabAlpha, dv, HT.data(), nNp);
    addFseep(dNp, nNp, ndm, c.kbar.data(), c.rhoF, c.drive.data(), dv, FS.data());
  }

  // dofMap vs stored DOFMAP
  int uOff[16], pOff[16];
  int total = dofMap(nNodes, ndm, pCarrier, uOff, pOff);
  bool dofOK = (total == c.totalDof);
  for (int n = 0; n < nNodes; n++) {
    if (uOff[n] != c.uOff[n] || pOff[n] != c.pOff[n]) dofOK = false;
  }
  if (!dofOK) g_fail++;
  printf("  [%s] DOFMAP total=%d(want %d)  %s\n", c.name.c_str(), total,
         c.totalDof, dofOK ? "ok" : "FAIL");

  compare("Q", c, Q, c.Q.v);
  compare("H", c, H, c.H.v);
  compare("S", c, S, c.S.v);
  compare("HT", c, HT, c.HT.v);
  compare("FSEEP", c, FS, c.FSEEP);
}

static void runCase(const Case& c) {
  printf("CASE %s (%s, %s)\n", c.name.c_str(), c.shape.c_str(), c.porder.c_str());
  if (c.shape == "T3")          driveCase<shapes::T3>(c);
  else if (c.shape == "Q4")     driveCase<shapes::Q4>(c);
  else if (c.shape == "H8")     driveCase<shapes::H8>(c);
  else if (c.shape == "BT6")    driveCase<shapes::BT6>(c);
  else if (c.shape == "BTET10") driveCase<shapes::BTET10>(c);
  else {
    printf("  UNKNOWN SHAPE %s\n", c.shape.c_str());
    g_fail++;
  }
}

// ---------------------------------------------------------------- cases file -- //
static int parseCases(const std::string& path) {
  std::ifstream f(path);
  if (!f) {
    printf("[cases] '%s' not present — self-tests only "
           "(sibling emits it; WP0.D cross-checks)\n", path.c_str());
    return 0;
  }
  int nCases = 0;
  std::string tok;
  while (f >> tok) {
    if (tok != "CASE") continue;
    Case c;
    std::string kw;
    f >> c.name >> kw >> c.shape >> kw >> c.porder >> kw >> c.ndm >> kw >> c.nNodes;
    while (f >> tok && tok != "END") {
      if (tok == "COORDS") {
        c.coords.resize((size_t)c.nNodes * c.ndm);
        for (auto& x : c.coords) f >> x;
      } else if (tok == "PARAMS") {
        f >> c.biotAlpha >> c.oneOverQbar >> c.rhoF >> c.stabAlpha >> c.thick;
        c.kbar.resize(c.ndm);
        for (auto& x : c.kbar) f >> x;
      } else if (tok == "DRIVE") {
        c.drive.resize(c.ndm);
        for (auto& x : c.drive) f >> x;
      } else if (tok == "DOFMAP") {
        f >> c.totalDof;
        c.uOff.resize(c.nNodes);
        c.pOff.resize(c.nNodes);
        for (auto& x : c.uOff) f >> x;
        for (auto& x : c.pOff) f >> x;
      } else if (tok == "Q") {
        c.Q = readMat(f);
      } else if (tok == "H") {
        c.H = readMat(f);
      } else if (tok == "S") {
        c.S = readMat(f);
      } else if (tok == "HT") {
        c.HT = readMat(f);
      } else if (tok == "FSEEP") {
        int len;
        f >> len;
        c.FSEEP.resize(len);
        for (auto& x : c.FSEEP) f >> x;
      }
    }
    runCase(c);
    nCases++;
  }
  printf("[cases] parsed %d case(s) from %s\n", nCases, path.c_str());
  return nCases;
}

// ----------------------------------------------------------------- self-tests // //
static void chk(const char* what, bool cond, double got, double want) {
  if (!cond) {
    printf("  SELFTEST FAIL %-28s got=%.15g want=%.15g\n", what, got, want);
    g_fail++;
  } else {
    printf("  selftest ok   %-28s (%.12g)\n", what, got);
  }
}

template <class P>
static void providerInvariants(const char* nm, const double* xy, bool hasTH,
                               double volume) {
  double volAcc = 0.0;
  bool puOK = true, gradOK = true, detOK = true, thOK = true;
  for (int gp = 0; gp < P::nGP; gp++) {
    double Nu[16], dNu[48], detJ = 0.0;
    P::eval(gp, xy, Nu, dNu, &detJ);
    if (detJ <= 0.0) detOK = false;
    volAcc += P::gpWeight[gp] * detJ;
    double sN = 0.0;
    double sG[3] = {0, 0, 0};
    for (int a = 0; a < P::nNodes; a++) {
      sN += Nu[a];
      for (int i = 0; i < P::ndm; i++) sG[i] += dNu[a * P::ndm + i];
    }
    if (std::fabs(sN - 1.0) > 1e-12) puOK = false;
    for (int i = 0; i < P::ndm; i++)
      if (std::fabs(sG[i]) > 1e-9) gradOK = false;
    // TH pressure partition (compacted linear barycentric)
    if (hasTH) {
      double Np[16], dNp[48];
      P::evalP(gp, xy, true, Np, dNp);
      int nNp = 0;
      for (int n = 0; n < P::nNodes; n++) nNp += P::pCarrierTH[n];
      double sP = 0.0, sPG[3] = {0, 0, 0};
      for (int b = 0; b < nNp; b++) {
        sP += Np[b];
        for (int i = 0; i < P::ndm; i++) sPG[i] += dNp[b * P::ndm + i];
      }
      if (std::fabs(sP - 1.0) > 1e-12) thOK = false;
      for (int i = 0; i < P::ndm; i++)
        if (std::fabs(sPG[i]) > 1e-9) thOK = false;
    }
  }
  char b[96];
  snprintf(b, sizeof b, "%s partition-of-unity", nm);   chk(b, puOK, 0, 0);
  snprintf(b, sizeof b, "%s Sum(dN)=0", nm);            chk(b, gradOK, 0, 0);
  snprintf(b, sizeof b, "%s detJ>0", nm);               chk(b, detOK, 0, 0);
  snprintf(b, sizeof b, "%s Sum(w*detJ)=vol", nm);
  chk(b, std::fabs(volAcc - volume) < 1e-10, volAcc, volume);
  if (hasTH) {
    snprintf(b, sizeof b, "%s TH-p partition", nm);     chk(b, thOK, 0, 0);
  }
}

static void selfTests() {
  printf("SELFTESTS\n");

  // --- provider invariants on positively-oriented canonical geometries ------ //
  double xyT3[6] = {0, 0,  1, 0,  0, 1};
  double xyQ4[8] = {0, 0,  1, 0,  1, 1,  0, 1};
  double xyH8[24] = {0, 0, 0,  1, 0, 0,  1, 1, 0,  0, 1, 0,
                     0, 0, 1,  1, 0, 1,  1, 1, 1,  0, 1, 1};
  // BT6: vertices (0,0),(1,0),(0,1); mid-edge nodes at midpoints (straight).
  double xyBT6[12] = {0, 0,  1, 0,  0, 1,   0.5, 0,  0.5, 0.5,  0, 0.5};
  // BTET10: positively-oriented tet vertices + straight-edge midpoints.
  double xyTet[30] = {
      0, 0, 0,   0, 1, 0,   1, 0, 0,   0, 0, 1,          // vertices 0..3
      0, 0.5, 0,   0.5, 0.5, 0,   0.5, 0, 0,             // 4:(0-1) 5:(1-2) 6:(0-2)
      0, 0, 0.5,   0.5, 0, 0.5,   0, 0.5, 0.5};          // 7:(0-3) 8:(2-3) 9:(1-3)

  providerInvariants<shapes::T3>("T3", xyT3, false, 0.5);
  providerInvariants<shapes::Q4>("Q4", xyQ4, false, 1.0);
  providerInvariants<shapes::H8>("H8", xyH8, false, 1.0);
  providerInvariants<shapes::BT6>("BT6", xyBT6, true, 0.5);
  providerInvariants<shapes::BTET10>("BTET10", xyTet, true, 1.0 / 6.0);

  // --- elemSizeH (largest vertex-pair distance) ----------------------------- //
  chk("T3 elemSizeH", std::fabs(shapes::T3::elemSizeH(xyT3) - std::sqrt(2.0)) < 1e-12,
      shapes::T3::elemSizeH(xyT3), std::sqrt(2.0));                 // hypotenuse
  chk("Q4 elemSizeH", std::fabs(shapes::Q4::elemSizeH(xyQ4) - std::sqrt(2.0)) < 1e-12,
      shapes::Q4::elemSizeH(xyQ4), std::sqrt(2.0));                 // diagonal
  chk("H8 elemSizeH", std::fabs(shapes::H8::elemSizeH(xyH8) - std::sqrt(3.0)) < 1e-12,
      shapes::H8::elemSizeH(xyH8), std::sqrt(3.0));                 // body diagonal
  chk("BT6 elemSizeH", std::fabs(shapes::BT6::elemSizeH(xyBT6) - std::sqrt(2.0)) < 1e-12,
      shapes::BT6::elemSizeH(xyBT6), std::sqrt(2.0));               // vertex-only
  chk("BTET10 elemSizeH",
      std::fabs(shapes::BTET10::elemSizeH(xyTet) - std::sqrt(2.0)) < 1e-12,
      shapes::BTET10::elemSizeH(xyTet), std::sqrt(2.0));            // vertex-only

  // --- scalar kernel helpers vs hand-computed contract formulas ------------- //
  // oneOverQbar: n/Kf (+ (alpha-n)/Ks when Ks>0)
  {
    double v = oneOverQbar(0.4, 2.2e6, 1.0, -1.0);       // Ks<=0 => infinite grains
    chk("oneOverQbar (Ks<=0)", std::fabs(v - 0.4 / 2.2e6) < 1e-20, v, 0.4 / 2.2e6);
    double w = oneOverQbar(0.4, 2.2e6, 0.79, 1.0e7);
    double want = 0.4 / 2.2e6 + (0.79 - 0.4) / 1.0e7;
    chk("oneOverQbar (finite Ks)", std::fabs(w - want) < 1e-18, w, want);
  }
  // stabAlphaAuto = alpha0*h^2/(Ks + 4Gs/3)
  {
    double a = stabAlphaAuto(2.0, 1.0e5, 6.0e4, 0.25);
    double want = 0.25 * 4.0 / (1.0e5 + 4.0 * 6.0e4 / 3.0);
    chk("stabAlphaAuto", std::fabs(a - want) < 1e-18, a, want);
  }

  // --- dofMap: heterogeneous-ndf (Taylor-Hood) placement -------------------- //
  {
    int pc[6] = {1, 1, 1, 0, 0, 0};   // BT6 TH carriers
    int uO[6], pO[6];
    int tot = dofMap(6, 2, pc, uO, pO);
    // node0 u@0 p@2 | node1 u@3 p@5 | node2 u@6 p@8 | node3 u@9 (no p)
    //   node4 u@11 | node5 u@13 -> total 15
    bool ok = (tot == 15) && uO[0] == 0 && pO[0] == 2 && uO[1] == 3 && pO[1] == 5 &&
              uO[2] == 6 && pO[2] == 8 && uO[3] == 9 && pO[3] == -1 && uO[4] == 11 &&
              pO[4] == -1 && uO[5] == 13 && pO[5] == -1;
    chk("dofMap TH (BT6)", ok, tot, 15);
  }
}

int main(int argc, char** argv) {
  std::string path = (argc > 1) ? argv[1] : "tests/ladruno_up_cases.txt";
  parseCases(path);
  selfTests();
  if (g_fail) {
    printf("\nRESULT: FAIL (%d check(s) failed)\n", g_fail);
    return 1;
  }
  printf("\nRESULT: PASS\n");
  return 0;
}
