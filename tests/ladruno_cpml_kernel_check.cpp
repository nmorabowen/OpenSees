// Standalone C++17 verification harness for LadrunoCPMLShapes.h and
// LadrunoCPMLKernel.h (Phase P0 of LadrunoCPML).
// Zero OpenSees dependencies.
//
// Modes:
//   1. `ladruno_cpml_kernel_check --self-check`
//      Runs comprehensive standalone C++ self-checks on all 9 topologies:
//        - Shape function partition of unity, zero gradient sum, linear
//          coordinate/gradient reproduction, positive Jacobian, exact volume.
//        - Standard right-handed TET10/BTET10 orientation support (signed
//          detJ < 0 with positive w_vol = w * |detJ|).
//        - Finite-difference tangent consistency: dR_int/du = K_eff,
//          dR_int/dv = C_eff, dR_int/da = M_eff, and Newmark combined
//          dR_int/d(du) = K_T_dyn across all 9 topologies, box & plane modes,
//          consistent & HRZ-lumped mass, multi-step history.
//        - Symmetry and positive-definiteness of K_eff, C_eff, M_eff.
//        - Static limit (dt <= 0.0 and setStaticMode(true)), interior GP
//          limit (xi_m = 0), classical unshifted PML limit (alpha_max = 0),
//          input validation (isValid()), and direct mutation-killing physics
//          invariants (kappa_max, alpha_max, damping d_m, corner cross terms,
//          and recursive-convolution memory sign).
//        - revertToLastCommit() (including bit-identical packState before vs
//          after revertToLastCommit), revertToStart(), and packState/unpackState
//          round-trip verification.
//   2. `ladruno_cpml_kernel_check [path/to/ladruno_cpml_cases.txt]`
//      Runs `--self-check` AND parses every CASE and STEP in the oracle file,
//      comparing R_int, K_eff, C_eff, M_eff, K_T_dyn, and P_ext against
//      CPMLElementKernel to <= 1e-12.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "LadrunoCPMLShapes.h"
#include "LadrunoCPMLKernel.h"

using namespace LadrunoCPML;

static constexpr double ORACLE_TOL = 1.0e-12;
static int g_failures = 0;
static double g_globalMaxAbsErr = 0.0;
static double g_globalMaxRelErr = 0.0;

// --------------------------------------------------------------------------- //
// Tokenizer (strips '#' comments to end of line)
// --------------------------------------------------------------------------- //
struct TokenStream {
  std::vector<std::string> tokens;
  size_t pos = 0;

  bool loadFile(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) return false;
    std::string line;
    while (std::getline(ifs, line)) {
      const size_t hashPos = line.find('#');
      if (hashPos != std::string::npos) {
        line.erase(hashPos);
      }
      std::istringstream iss(line);
      std::string tok;
      while (iss >> tok) {
        tokens.push_back(tok);
      }
    }
    pos = 0;
    return true;
  }

  bool hasMore() const { return pos < tokens.size(); }

  bool peekStr(std::string& out) const {
    if (pos >= tokens.size()) return false;
    out = tokens[pos];
    return true;
  }

  bool nextStr(std::string& out) {
    if (pos >= tokens.size()) return false;
    out = tokens[pos++];
    return true;
  }

  bool expect(const char* kw) {
    std::string s;
    if (!nextStr(s) || s != kw) {
      std::printf("PARSE ERROR: expected '%s', got '%s' at token %zu\n",
                  kw, s.c_str(), pos);
      return false;
    }
    return true;
  }

  bool nextInt(int& out) {
    std::string s;
    if (!nextStr(s)) return false;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str() || *endp != '\0') return false;
    out = static_cast<int>(v);
    return true;
  }

  bool nextDouble(double& out) {
    std::string s;
    if (!nextStr(s)) return false;
    char* endp = nullptr;
    out = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str() || *endp != '\0') return false;
    return true;
  }

  bool nextVec(std::vector<double>& v, size_t count) {
    v.resize(count);
    for (size_t i = 0; i < count; ++i) {
      if (!nextDouble(v[i])) return false;
    }
    return true;
  }
};

static bool parseTopologyName(const std::string& s, CPMLTopology& topo) {
  if (s == "T3")          { topo = CPMLTopology::T3; return true; }
  if (s == "Q4")          { topo = CPMLTopology::Q4; return true; }
  if (s == "T6")          { topo = CPMLTopology::T6; return true; }
  if (s == "BT6")         { topo = CPMLTopology::BT6; return true; }
  if (s == "TET4")        { topo = CPMLTopology::TET4; return true; }
  if (s == "H8")          { topo = CPMLTopology::H8; return true; }
  if (s == "TET10")       { topo = CPMLTopology::TET10; return true; }
  if (s == "BTET10")      { topo = CPMLTopology::BTET10; return true; }
  if (s == "H20")         { topo = CPMLTopology::H20; return true; }
  return false;
}

// --------------------------------------------------------------------------- //
// Comparison helper
// --------------------------------------------------------------------------- //
static void compareArray(const std::string& caseName, int stepIdx,
                         const char* label,
                         const std::vector<double>& got,
                         const std::vector<double>& want,
                         double& caseMaxAbs, double& caseMaxRel,
                         bool& caseOK) {
  if (got.size() != want.size()) {
    std::printf("  [%s step %d] %s SIZE MISMATCH got=%zu want=%zu FAIL\n",
                caseName.c_str(), stepIdx, label, got.size(), want.size());
    caseOK = false;
    g_failures++;
    return;
  }
  double maxWant = 0.0;
  for (double w : want) {
    if (std::isfinite(w)) maxWant = std::max(maxWant, std::fabs(w));
  }
  const double scale = std::max(1.0, maxWant);

  double maxAbs = 0.0;
  bool localOK = true;
  for (size_t i = 0; i < got.size(); ++i) {
    if (!std::isfinite(got[i]) || !std::isfinite(want[i])) {
      localOK = false;
      continue;
    }
    const double diff = std::fabs(got[i] - want[i]);
    if (diff > maxAbs) maxAbs = diff;
    if (diff / scale > ORACLE_TOL) {
      localOK = false;
    }
  }
  const double maxRel = maxAbs / scale;
  caseMaxAbs = std::max(caseMaxAbs, maxAbs);
  caseMaxRel = std::max(caseMaxRel, maxRel);
  g_globalMaxAbsErr = std::max(g_globalMaxAbsErr, maxAbs);
  g_globalMaxRelErr = std::max(g_globalMaxRelErr, maxRel);

  if (!localOK) {
    std::printf("  [%s step %d] %-8s FAIL maxAbs=%.3e maxRel=%.3e (scale=%.3e)\n",
                caseName.c_str(), stepIdx, label, maxAbs, maxRel, scale);
    caseOK = false;
    g_failures++;
  }
}

static bool isPressureDOF(const CPMLElementKernel& kernel, int dof) {
  const int nP = kernel.numPressureDOFs();
  for (int b = 0; b < nP; ++b) {
    if (kernel.pOffset(b) == dof) return true;
  }
  return false;
}

static void compareUPVec(const std::string& caseName, int stepIdx,
                         const char* label, const CPMLElementKernel& kernel,
                         const std::vector<double>& got,
                         const std::vector<double>& want,
                         double& caseMaxAbs, double& caseMaxRel,
                         bool& caseOK) {
  compareArray(caseName, stepIdx, label, got, want, caseMaxAbs, caseMaxRel, caseOK);
  if (!kernel.isUPMode() || got.size() != want.size()) return;
  const int nDOF = kernel.numDOF();
  double scaleU = 0.0, scaleP = 0.0;
  for (int i = 0; i < nDOF; ++i) {
    if (isPressureDOF(kernel, i)) {
      scaleP = std::max(scaleP, std::fabs(want[i]));
    } else {
      scaleU = std::max(scaleU, std::fabs(want[i]));
    }
  }
  if (scaleU <= 0.0) scaleU = 1.0;
  if (scaleP <= 0.0) scaleP = 1.0;
  for (int i = 0; i < nDOF; ++i) {
    const double sc = isPressureDOF(kernel, i) ? scaleP : scaleU;
    const double rel = std::fabs(got[i] - want[i]) / sc;
    caseMaxRel = std::max(caseMaxRel, rel);
    g_globalMaxRelErr = std::max(g_globalMaxRelErr, rel);
    if (rel > 1.0e-11) {
      std::printf("  [%s step %d] %-8s UP-BLOCK FAIL dof=%d rel=%.3e (sc=%.3e)\n",
                  caseName.c_str(), stepIdx, label, i, rel, sc);
      caseOK = false;
      g_failures++;
      break;
    }
  }
}

static void compareUPMat(const std::string& caseName, int stepIdx,
                         const char* label, const CPMLElementKernel& kernel,
                         const std::vector<double>& got,
                         const std::vector<double>& want,
                         double& caseMaxAbs, double& caseMaxRel,
                         bool& caseOK) {
  compareArray(caseName, stepIdx, label, got, want, caseMaxAbs, caseMaxRel, caseOK);
  if (!kernel.isUPMode() || got.size() != want.size()) return;
  const int nDOF = kernel.numDOF();
  double sc[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  for (int i = 0; i < nDOF; ++i) {
    const int bi = isPressureDOF(kernel, i) ? 1 : 0;
    for (int j = 0; j < nDOF; ++j) {
      const int bj = isPressureDOF(kernel, j) ? 1 : 0;
      sc[bi][bj] = std::max(sc[bi][bj], std::fabs(want[(size_t)i * nDOF + j]));
    }
  }
  for (int bi = 0; bi < 2; ++bi) {
    for (int bj = 0; bj < 2; ++bj) {
      if (sc[bi][bj] <= 0.0) sc[bi][bj] = 1.0;
    }
  }
  for (int i = 0; i < nDOF; ++i) {
    const int bi = isPressureDOF(kernel, i) ? 1 : 0;
    for (int j = 0; j < nDOF; ++j) {
      const int bj = isPressureDOF(kernel, j) ? 1 : 0;
      const double diff = std::fabs(got[(size_t)i * nDOF + j] - want[(size_t)i * nDOF + j]);
      const double rel = diff / sc[bi][bj];
      caseMaxRel = std::max(caseMaxRel, rel);
      g_globalMaxRelErr = std::max(g_globalMaxRelErr, rel);
      if (rel > 1.0e-11) {
        std::printf("  [%s step %d] %-8s UP-BLOCK(%d,%d) FAIL (%d,%d) rel=%.3e (sc=%.3e)\n",
                    caseName.c_str(), stepIdx, label, bi, bj, i, j, rel, sc[bi][bj]);
        caseOK = false;
        g_failures++;
        return;
      }
    }
  }
}

// --------------------------------------------------------------------------- //
// Oracle Cases Runner
// --------------------------------------------------------------------------- //
static int runCasesFile(const std::string& path) {
  TokenStream ts;
  if (!ts.loadFile(path)) {
    std::printf("[cases] Could not open '%s'\n", path.c_str());
    g_failures++;
    return 0;
  }

  int nCases = 0;
  while (ts.hasMore()) {
    std::string tok;
    if (!ts.nextStr(tok)) break;
    if (tok != "CASE") {
      std::printf("PARSE ERROR: expected 'CASE', got '%s'\n", tok.c_str());
      g_failures++;
      break;
    }

    std::string caseName, topoStr, geomModeStr;
    int ndm = 0, nNodes = 0, lumpedInt = 0, nGeomParams = 0, nSteps = 0;
    double E = 0.0, nu = 0.0, rho = 0.0;
    double m_deg = 0.0, R0 = 0.0, kappa_max = 0.0, alpha_max = 0.0;
    double c1 = 0.0, c2 = 0.0, c3 = 0.0;
    std::vector<double> geomParams, coords;
    CPMLBiotParams biot{};

    bool ok = ts.nextStr(caseName) &&
              ts.expect("Topology") && ts.nextStr(topoStr) &&
              ts.nextInt(ndm) && ts.nextInt(nNodes) && ts.nextInt(lumpedInt) &&
              ts.expect("Material") && ts.nextDouble(E) && ts.nextDouble(nu) && ts.nextDouble(rho);

    std::string nextKw;
    if (ok && ts.peekStr(nextKw) && nextKw == "BiotParams") {
      ts.nextStr(nextKw);
      biot.enabled = true;
      ok = ts.nextDouble(biot.Kf) && ts.nextDouble(biot.poro) &&
           ts.nextDouble(biot.rhoF) && ts.nextDouble(biot.biotAlpha) &&
           ts.nextDouble(biot.Ks) && ts.nextInt(biot.pOrder) &&
           ts.nextInt(biot.formulation) && ts.nextInt(biot.stabMode) &&
           ts.nextDouble(biot.stabValue);
      for (int m = 0; ok && m < ndm && m < 3; ++m) {
        ok = ts.nextDouble(biot.kbar[m]);
      }
    }

    ok = ok &&
         ts.expect("PMLParams") && ts.nextDouble(m_deg) && ts.nextDouble(R0) &&
         ts.nextDouble(kappa_max) && ts.nextDouble(alpha_max) &&
         ts.expect("GeomMode") && ts.nextStr(geomModeStr) && ts.nextInt(nGeomParams) &&
         ts.nextVec(geomParams, (size_t)nGeomParams) &&
         ts.expect("Coords") && ts.nextVec(coords, (size_t)(nNodes * ndm)) &&
         ts.expect("NewmarkCoeffs") && ts.nextDouble(c1) && ts.nextDouble(c2) && ts.nextDouble(c3) &&
         ts.expect("NumSteps") && ts.nextInt(nSteps);

    if (!ok) {
      std::printf("CASE %s: header parse failed\n", caseName.c_str());
      g_failures++;
      break;
    }

    CPMLTopology topo;
    if (!parseTopologyName(topoStr, topo)) {
      std::printf("CASE %s: unknown topology '%s'\n", caseName.c_str(), topoStr.c_str());
      g_failures++;
      break;
    }

    CPMLGeomConfig geom;
    if (geomModeStr == "box") {
      if (ndm == 2 && nGeomParams == 6) {
        geom = CPMLGeomConfig::makeBox2D(geomParams[0], geomParams[1],
                                         geomParams[2], geomParams[3],
                                         geomParams[4], geomParams[5]);
      } else if (ndm == 3 && nGeomParams == 9) {
        geom = CPMLGeomConfig::makeBox3D(geomParams[0], geomParams[1],
                                         geomParams[2], geomParams[3],
                                         geomParams[4], geomParams[5],
                                         geomParams[6], geomParams[7],
                                         geomParams[8]);
      } else {
        std::printf("CASE %s: invalid box geom params count %d for ndm=%d\n",
                    caseName.c_str(), nGeomParams, ndm);
        g_failures++;
        break;
      }
    } else if (geomModeStr == "plane") {
      if (ndm == 2 && nGeomParams == 5) {
        geom = CPMLGeomConfig::makePlane2D(geomParams[0], geomParams[1],
                                           geomParams[2], geomParams[3],
                                           geomParams[4]);
      } else if (ndm == 3 && nGeomParams == 7) {
        geom = CPMLGeomConfig::makePlane3D(geomParams[0], geomParams[1], geomParams[2],
                                           geomParams[3], geomParams[4], geomParams[5],
                                           geomParams[6]);
      } else {
        std::printf("CASE %s: invalid plane geom params count %d for ndm=%d\n",
                    caseName.c_str(), nGeomParams, ndm);
        g_failures++;
        break;
      }
    } else {
      std::printf("CASE %s: unknown GeomMode '%s'\n", caseName.c_str(), geomModeStr.c_str());
      g_failures++;
      break;
    }

    CPMLMaterialParams mat(E, nu, rho);
    CPMLGradingParams grading(m_deg, R0, kappa_max, alpha_max);
    CPMLElementKernel kernel;
    if (!kernel.init(topo, lumpedInt != 0, coords.data(), mat, grading, geom, biot)) {
      std::printf("CASE %s: kernel.init reported invalid geometry/parameters\n", caseName.c_str());
      g_failures++;
    }

    const size_t nDOF = (size_t)kernel.numDOF();
    const size_t nMat = nDOF * nDOF;
    double caseMaxAbs = 0.0, caseMaxRel = 0.0;
    bool caseOK = true;

    for (int s = 0; s < nSteps; ++s) {
      int stepIdx = -1;
      double dt = 0.0;
      std::vector<double> u, v, a, wantR, wantK, wantC, wantM, wantKT, wantP(1);
      bool stepParseOK =
          ts.expect("STEP") && ts.nextInt(stepIdx) && ts.nextDouble(dt) &&
          ts.expect("U") && ts.nextVec(u, nDOF) &&
          ts.expect("V") && ts.nextVec(v, nDOF) &&
          ts.expect("A") && ts.nextVec(a, nDOF) &&
          ts.expect("R_int") && ts.nextVec(wantR, nDOF) &&
          ts.expect("K_eff") && ts.nextVec(wantK, nMat) &&
          ts.expect("C_eff") && ts.nextVec(wantC, nMat) &&
          ts.expect("M_eff") && ts.nextVec(wantM, nMat) &&
          ts.expect("K_T_dyn") && ts.nextVec(wantKT, nMat) &&
          ts.expect("P_ext") && ts.nextDouble(wantP[0]) &&
          ts.expect("END_STEP");

      if (!stepParseOK || stepIdx != s) {
        std::printf("CASE %s: step %d parse failed\n", caseName.c_str(), s);
        g_failures++;
        caseOK = false;
        break;
      }

      kernel.setTrial(u.data(), v.data(), a.data(), dt);

      std::vector<double> gotP = {kernel.getP_ext()};
      std::vector<double> gotKT = kernel.getK_T_dyn(c1, c2, c3);

      compareUPVec(caseName, s, "R_int",   kernel, kernel.getR_int(), wantR,  caseMaxAbs, caseMaxRel, caseOK);
      compareUPMat(caseName, s, "K_eff",   kernel, kernel.getK_eff(), wantK,  caseMaxAbs, caseMaxRel, caseOK);
      compareUPMat(caseName, s, "C_eff",   kernel, kernel.getC_eff(), wantC,  caseMaxAbs, caseMaxRel, caseOK);
      compareUPMat(caseName, s, "M_eff",   kernel, kernel.getM_eff(), wantM,  caseMaxAbs, caseMaxRel, caseOK);
      compareUPMat(caseName, s, "K_T_dyn", kernel, gotKT,             wantKT, caseMaxAbs, caseMaxRel, caseOK);
      compareArray(caseName, s, "P_ext",   gotP,              wantP,  caseMaxAbs, caseMaxRel, caseOK);

      kernel.commit();
    }

    if (!ts.expect("END_CASE")) {
      g_failures++;
      break;
    }

    std::printf("CASE %-28s steps=%d maxAbs=%.3e maxRel=%.3e %s\n",
                caseName.c_str(), nSteps, caseMaxAbs, caseMaxRel,
                caseOK ? "PASS" : "FAIL");
    nCases++;
  }

  std::printf("[cases] verified %d case(s) from %s (global maxAbs=%.3e, maxRel=%.3e)\n",
              nCases, path.c_str(), g_globalMaxAbsErr, g_globalMaxRelErr);
  return nCases;
}

// --------------------------------------------------------------------------- //
// Standalone Self-Checks
// --------------------------------------------------------------------------- //
static void recordCheck(const char* name, bool ok, double err = 0.0) {
  if (!ok) {
    std::printf("  [SELF-CHECK FAIL] %-46s (err=%.3e)\n", name, err);
    g_failures++;
  } else {
    std::printf("  [self-check ok]   %-46s (err=%.3e)\n", name, err);
  }
}

static void checkTopologyShapes(const char* name, CPMLTopology topo,
                                const std::vector<double>& coords,
                                double expectedVol,
                                bool allowNegDetJ = false) {
  const int ndm = shapes::topologyNdm(topo);
  const int nNodes = shapes::topologyNumNodes(topo);
  const int nGP = shapes::topologyNumGP(topo);

  double vol = 0.0;
  double maxPuErr = 0.0;
  double maxGradSumErr = 0.0;
  double maxLinGradErr = 0.0;
  int nPos = 0, nNeg = 0;

  for (int g = 0; g < nGP; ++g) {
    double N[shapes::MAX_NODES]{};
    double dNdx[shapes::MAX_NODES * 3]{};
    double x_gp[3]{};
    double detJ = 0.0, w_vol = 0.0;
    shapes::evalTopologyGP(topo, g, coords.data(), N, dNdx, &detJ, x_gp, &w_vol);
    if (detJ > 1.0e-18) nPos++;
    else if (detJ < -1.0e-18) nNeg++;
    vol += w_vol;

    double sumN = 0.0;
    double sumGrad[3] = {0.0, 0.0, 0.0};
    double gradX[3][3]{}; // gradX[i][j] = sum_a x_{a,i} dN_a/dx_j = delta_{ij}
    for (int a = 0; a < nNodes; ++a) {
      sumN += N[a];
      for (int j = 0; j < ndm; ++j) {
        sumGrad[j] += dNdx[a * ndm + j];
        for (int i = 0; i < ndm; ++i) {
          gradX[i][j] += coords[a * ndm + i] * dNdx[a * ndm + j];
        }
      }
    }
    maxPuErr = std::max(maxPuErr, std::fabs(sumN - 1.0));
    for (int j = 0; j < ndm; ++j) {
      maxGradSumErr = std::max(maxGradSumErr, std::fabs(sumGrad[j]));
      for (int i = 0; i < ndm; ++i) {
        const double want = (i == j) ? 1.0 : 0.0;
        maxLinGradErr = std::max(maxLinGradErr, std::fabs(gradX[i][j] - want));
      }
    }
  }

  const bool detOK = allowNegDetJ ? (nPos == nGP || nNeg == nGP) : (nPos == nGP);
  const double volErr = std::fabs(vol - expectedVol);
  const double worstErr = std::max({maxPuErr, maxGradSumErr, maxLinGradErr, volErr});
  char label[96];
  std::snprintf(label, sizeof(label), "Shapes %s (PU, grad, detJ, vol)", name);
  recordCheck(label, detOK && worstErr < 1.0e-12, worstErr);
}

static bool isSymmetricAndSPD(const std::vector<double>& M, int n,
                              double& symErr, double& minQuad) {
  double maxAbsEntry = 0.0;
  for (double v : M) maxAbsEntry = std::max(maxAbsEntry, std::fabs(v));
  const double scale = std::max(1.0, maxAbsEntry);
  symErr = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      symErr = std::max(symErr, std::fabs(M[(size_t)i * n + j] - M[(size_t)j * n + i]) / scale);
    }
  }
  minQuad = 1.0e300;
  for (int i = 0; i < n; ++i) {
    minQuad = std::min(minQuad, M[(size_t)i * n + i]);
  }
  std::vector<double> z(n);
  for (int trial = 1; trial <= 4; ++trial) {
    for (int i = 0; i < n; ++i) {
      z[i] = std::sin(0.73 * (i + 1) * trial + 0.19 * trial);
    }
    double q = 0.0;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        q += z[i] * M[(size_t)i * n + j] * z[j];
    minQuad = std::min(minQuad, q);
  }
  return (symErr < 1.0e-14) && (minQuad > 0.0);
}

static void checkKernelFiniteDiffAndState(const char* label,
                                          CPMLTopology topo,
                                          bool lumped,
                                          const std::vector<double>& coords,
                                          const CPMLGeomConfig& geom) {
  const CPMLMaterialParams mat(2.5e7, 0.25, 2000.0);
  const CPMLGradingParams grading(2.0, 1.0e-6, 3.0, 15.0);

  CPMLElementKernel kernel(topo, lumped, coords.data(), mat, grading, geom);
  if (!kernel.isValid()) {
    recordCheck(label, false, 1.0);
    return;
  }
  const int nDOF = kernel.numDOF();

  // Step 0: advance one step and commit so recursive convolution memory is active
  const double dt = 0.005;
  std::vector<double> u0(nDOF), v0(nDOF), a0(nDOF);
  for (int i = 0; i < nDOF; ++i) {
    u0[i] = 1.2e-4 * std::sin(0.37 * (i + 1));
    v0[i] = 2.5e-2 * std::cos(0.51 * (i + 1));
    a0[i] = -1.8e0 * std::sin(0.29 * (i + 1));
  }
  kernel.setTrial(u0.data(), v0.data(), a0.data(), dt);
  kernel.commit();

  // Pack committed state immediately after commit()
  std::vector<double> packedAfterCommit(kernel.numStateDoubles(), 0.0);
  kernel.packState(packedAfterCommit.data());

  // Step 1: evaluate at a second non-trivial trial state
  std::vector<double> u1(nDOF), v1(nDOF), a1(nDOF);
  for (int i = 0; i < nDOF; ++i) {
    u1[i] = u0[i] + 0.8e-4 * std::cos(0.43 * (i + 1));
    v1[i] = v0[i] - 1.4e-2 * std::sin(0.61 * (i + 1));
    a1[i] = a0[i] + 2.1e0  * std::cos(0.17 * (i + 1));
  }
  kernel.setTrial(u1.data(), v1.data(), a1.data(), dt);

  const double c1 = 1.0;
  const double c2 = 2.0 / dt;
  const double c3 = 4.0 / (dt * dt);
  const std::vector<double> R_base = kernel.getR_int();
  const std::vector<double> K_eff  = kernel.getK_eff();
  const std::vector<double> C_eff  = kernel.getC_eff();
  const std::vector<double> M_eff  = kernel.getM_eff();
  const std::vector<double> KT_dyn = kernel.getK_T_dyn(c1, c2, c3);

  // Check SPD of K_eff, C_eff, M_eff
  double symK = 0.0, minK = 0.0, symC = 0.0, minC = 0.0, symM = 0.0, minM = 0.0;
  const bool spdOK = isSymmetricAndSPD(K_eff, nDOF, symK, minK) &&
                     isSymmetricAndSPD(C_eff, nDOF, symC, minC) &&
                     isSymmetricAndSPD(M_eff, nDOF, symM, minM);
  char spdMsg[96];
  std::snprintf(spdMsg, sizeof(spdMsg), "SPD K/C/M %s", label);
  recordCheck(spdMsg, spdOK, std::max({symK, symC, symM}));

  // Central finite-difference check of K_T_dyn = c1*dR/du + c2*dR/dv + c3*dR/da
  const double eps = 1.0e-7;
  double maxFdRelErr = 0.0;
  double maxKT = 0.0;
  for (double val : KT_dyn) maxKT = std::max(maxKT, std::fabs(val));

  for (int j = 0; j < nDOF; ++j) {
    std::vector<double> up = u1, vp = v1, ap = a1;
    std::vector<double> um = u1, vm = v1, am = a1;
    up[j] += c1 * eps; vp[j] += c2 * eps; ap[j] += c3 * eps;
    um[j] -= c1 * eps; vm[j] -= c2 * eps; am[j] -= c3 * eps;

    kernel.setTrial(up.data(), vp.data(), ap.data(), dt);
    const std::vector<double> Rp = kernel.getR_int();

    kernel.setTrial(um.data(), vm.data(), am.data(), dt);
    const std::vector<double> Rm = kernel.getR_int();

    for (int i = 0; i < nDOF; ++i) {
      const double fd_ij = (Rp[i] - Rm[i]) / (2.0 * eps);
      const double kt_ij = KT_dyn[(size_t)i * nDOF + j];
      const double rel = std::fabs(fd_ij - kt_ij) / std::max(1.0, maxKT);
      maxFdRelErr = std::max(maxFdRelErr, rel);
    }
  }
  char fdMsg[96];
  std::snprintf(fdMsg, sizeof(fdMsg), "FD K_T_dyn %s", label);
  recordCheck(fdMsg, maxFdRelErr < 1.0e-6, maxFdRelErr);

  // Check revertToLastCommit():
  // (a) packState() immediately after revertToLastCommit() must be bit-identical
  //     to packedAfterCommit!
  kernel.revertToLastCommit();
  std::vector<double> packedAfterRevert(kernel.numStateDoubles(), 0.0);
  kernel.packState(packedAfterRevert.data());
  double revertPackDiff = 0.0;
  for (size_t i = 0; i < packedAfterCommit.size(); ++i) {
    revertPackDiff = std::max(revertPackDiff,
                              std::fabs(packedAfterRevert[i] - packedAfterCommit[i]));
  }

  // (b) Re-evaluating Step 1 after revertToLastCommit() reproduces R_base identically
  kernel.setTrial(u1.data(), v1.data(), a1.data(), dt);
  const std::vector<double> R_re = kernel.getR_int();
  double revErr = 0.0;
  for (int i = 0; i < nDOF; ++i) {
    revErr = std::max(revErr, std::fabs(R_re[i] - R_base[i]));
  }

  // Check unpackState into a fresh kernel reproduces Step 1 to 0.0 error
  CPMLElementKernel clone(topo, lumped, coords.data(), mat, grading, geom);
  clone.unpackState(packedAfterCommit.data());
  clone.setTrial(u1.data(), v1.data(), a1.data(), dt);
  const std::vector<double> R_clone = clone.getR_int();
  double packErr = 0.0;
  for (int i = 0; i < nDOF; ++i) {
    packErr = std::max(packErr, std::fabs(R_clone[i] - R_base[i]));
  }

  // Check revertToStart() zeroes memory and matches fresh kernel at Step 0
  kernel.revertToStart();
  kernel.setTrial(u0.data(), v0.data(), a0.data(), dt);
  CPMLElementKernel fresh(topo, lumped, coords.data(), mat, grading, geom);
  fresh.setTrial(u0.data(), v0.data(), a0.data(), dt);
  double startErr = 0.0;
  for (int i = 0; i < nDOF; ++i) {
    startErr = std::max(startErr, std::fabs(kernel.getR_int()[i] - fresh.getR_int()[i]));
  }

  char stateMsg[96];
  std::snprintf(stateMsg, sizeof(stateMsg), "State/Revert/Pack %s", label);
  const double worstStateErr = std::max({revertPackDiff, revErr, packErr, startErr});
  recordCheck(stateMsg, worstStateErr == 0.0, worstStateErr);
}

static void checkPhysicsLimitsAndMutations(const std::vector<double>& xyQ4) {
  const CPMLMaterialParams mat(2.5e7, 0.25, 2000.0);
  const CPMLGradingParams grading(2.0, 1.0e-6, 3.0, 15.0);
  const CPMLGeomConfig cornerBox = CPMLGeomConfig::makeBox2D(-1.0, 0.05, -1.0, 0.08, 1.5, 1.2);

  // 1. Static limit (dt = 0.0 and setStaticMode(true))
  {
    CPMLElementKernel k(CPMLTopology::Q4, false, xyQ4.data(), mat, grading, cornerBox);
    const int nDOF = k.numDOF();
    std::vector<double> u(nDOF), v(nDOF), a(nDOF);
    for (int i = 0; i < nDOF; ++i) {
      u[i] = 1.0e-3 * (i + 1);
      v[i] = 0.2 * (i + 1);
      a[i] = 5.0 * (i + 1);
    }
    // Advance one dynamic step and commit
    k.setTrial(u.data(), v.data(), a.data(), 0.005);
    k.commit();
    std::vector<double> stateBeforeStatic(k.numStateDoubles(), 0.0);
    k.packState(stateBeforeStatic.data());

    // Evaluate with dt = 0.0 (static)
    k.setTrial(u.data(), v.data(), a.data(), 0.0);
    const auto& K0 = k.getK0();
    const auto& Keff = k.getK_eff();
    const auto& Ceff = k.getC_eff();
    const auto& Rint = k.getR_int();

    double maxStaticErr = 0.0;
    for (size_t idx = 0; idx < K0.size(); ++idx) {
      maxStaticErr = std::max(maxStaticErr, std::fabs(Keff[idx] - K0[idx]));
      maxStaticErr = std::max(maxStaticErr, std::fabs(Ceff[idx]));
    }
    for (int i = 0; i < nDOF; ++i) {
      double wantRi = 0.0;
      for (int j = 0; j < nDOF; ++j) wantRi += K0[(size_t)i * nDOF + j] * u[j];
      maxStaticErr = std::max(maxStaticErr, std::fabs(Rint[i] - wantRi));
    }

    // Also test explicit setStaticMode(true) even when dt > 0
    k.setStaticMode(true);
    k.setTrial(u.data(), v.data(), a.data(), 0.005);
    for (int i = 0; i < nDOF; ++i) {
      double wantRi = 0.0;
      for (int j = 0; j < nDOF; ++j) wantRi += K0[(size_t)i * nDOF + j] * u[j];
      maxStaticErr = std::max(maxStaticErr, std::fabs(k.getR_int()[i] - wantRi));
    }
    k.setStaticMode(false);
    k.revertToLastCommit();
    std::vector<double> stateAfterStatic(k.numStateDoubles(), 0.0);
    k.packState(stateAfterStatic.data());
    for (size_t idx = 0; idx < stateBeforeStatic.size(); ++idx) {
      maxStaticErr = std::max(maxStaticErr,
                              std::fabs(stateAfterStatic[idx] - stateBeforeStatic[idx]));
    }
    recordCheck("Static limit (dt<=0 & staticMode -> K0*u)", maxStaticErr == 0.0, maxStaticErr);
  }

  // 2. Interior GP limit (xi_m = 0 -> K_eff == K_0, M_eff == M_0, C_eff == 0)
  //    and Classical PML limit (alpha_max = 0, d_m > 0)
  {
    const CPMLGeomConfig interiorBox = CPMLGeomConfig::makeBox2D(-10.0, 10.0, -10.0, 10.0, 1.5, 1.2);
    CPMLElementKernel kInt(CPMLTopology::Q4, false, xyQ4.data(), mat, grading, interiorBox);
    const int nDOF = kInt.numDOF();
    std::vector<double> u(nDOF, 1.0e-3), v(nDOF, 0.1), a(nDOF, 0.0);
    kInt.setTrial(u.data(), v.data(), a.data(), 0.005);

    double intErr = 0.0;
    for (size_t idx = 0; idx < kInt.getK0().size(); ++idx) {
      intErr = std::max(intErr, std::fabs(kInt.getK_eff()[idx] - kInt.getK0()[idx]));
      intErr = std::max(intErr, std::fabs(kInt.getM_eff()[idx] - kInt.getM0()[idx]));
      intErr = std::max(intErr, std::fabs(kInt.getC_eff()[idx]));
    }
    recordCheck("Interior GP limit (xi=0 -> K0, M0, C=0)", intErr < 1.0e-12, intErr);

    // Classical unshifted PML (alpha_max = 0.0, kappa_max = 1.0, d_m > 0)
    const CPMLGradingParams classicalGrading(2.0, 1.0e-6, 1.0, 0.0);
    CPMLElementKernel kClass(CPMLTopology::Q4, false, xyQ4.data(), mat, classicalGrading, cornerBox);
    kClass.setTrial(u.data(), v.data(), a.data(), 0.005);
    double symC = 0.0, minC = 0.0;
    const bool classOK = isSymmetricAndSPD(kClass.getC_eff(), nDOF, symC, minC);
    recordCheck("Classical unshifted PML limit (alpha=0 SPD)", classOK, symC);
  }

  // 3. Input validation (isValid() rejects bad material/grading/geom/inverted element)
  {
    const CPMLMaterialParams badMat1(-100.0, 0.25, 2000.0);
    const CPMLMaterialParams badMat2(2.5e7, 0.5, 2000.0);
    const CPMLGradingParams badGrad1(2.0, 1.2, 1.0, 0.0);
    const CPMLGradingParams badGrad2(2.0, 1.0e-6, 0.5, 0.0);
    const CPMLGeomConfig badGeom1 = CPMLGeomConfig::makePlane2D(0.0, 0.0, 0.0, 0.0, 1.0);
    const CPMLGeomConfig badGeom2 = CPMLGeomConfig::makeBox2D(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
    const std::vector<double> invertedQ4 = {0.0, 0.0,  0.0, 1.0,  1.0, 1.0,  1.0, 0.0};

    CPMLElementKernel kBad(CPMLTopology::Q4, false, invertedQ4.data(), mat, grading, cornerBox);
    const bool validChecksOK = !badMat1.isValid() && !badMat2.isValid() &&
                               !badGrad1.isValid() && !badGrad2.isValid() &&
                               !badGeom1.isValid() && !badGeom2.isValid() &&
                               !kBad.isValid();
    recordCheck("Input parameter & Jacobian validation", validChecksOK, 0.0);
  }

  // 4. Direct physics assertions that kill all 5 mutation modes:
  //    - CPML_NO_KAPPA: M_eff > M_0 when kappa_max > 1
  //    - CPML_NO_ALPHA: K_eff(alpha_max > 0) != K_eff(alpha_max = 0)
  //    - CPML_NO_DAMPING: C_eff > 0 when d_m > 0
  //    - CPML_NO_CORNER_CROSS: Rigid-body translation mode 1^T K_eff 1 > 0 in 2D corner
  //    - CPML_WRONG_SIGN_MEMORY: Exact 2-step RCFilter1D memory sign
  {
    CPMLElementKernel kFull(CPMLTopology::Q4, false, xyQ4.data(), mat, grading, cornerBox);
    const int nDOF = kFull.numDOF();
    std::vector<double> u(nDOF, 1.0e-3), v(nDOF, 0.1), a(nDOF, 0.0);
    kFull.setTrial(u.data(), v.data(), a.data(), 0.005);

    // (a) Kill CPML_NO_KAPPA: trace(M_eff) > trace(M_0)
    double trMeff = 0.0, trM0 = 0.0;
    for (int i = 0; i < nDOF; ++i) {
      trMeff += kFull.getM_eff()[(size_t)i * nDOF + i];
      trM0   += kFull.getM0()[(size_t)i * nDOF + i];
    }
    recordCheck("Physics: kappa_max>1 scales M_eff > M_0", trMeff > 1.5 * trM0, trMeff - trM0);

    // (b) Kill CPML_NO_DAMPING: trace(C_eff) > 0
    double trCeff = 0.0;
    for (int i = 0; i < nDOF; ++i) {
      trCeff += kFull.getC_eff()[(size_t)i * nDOF + i];
    }
    recordCheck("Physics: damping d_m>0 produces C_eff > 0", trCeff > 1.0e3, trCeff);

    // (c) Kill CPML_NO_ALPHA: compare K_eff with alpha_max=500 vs alpha_max=0
    const CPMLGradingParams gradAlpha0(2.0, 1.0e-6, 3.0, 0.0);
    const CPMLGradingParams gradAlphaHi(2.0, 1.0e-6, 3.0, 500.0);
    CPMLElementKernel kA0(CPMLTopology::Q4, false, xyQ4.data(), mat, gradAlpha0, cornerBox);
    CPMLElementKernel kAHi(CPMLTopology::Q4, false, xyQ4.data(), mat, gradAlphaHi, cornerBox);
    kA0.setTrial(u.data(), v.data(), a.data(), 0.005);
    kAHi.setTrial(u.data(), v.data(), a.data(), 0.005);
    double maxAlphaDiff = 0.0;
    for (size_t idx = 0; idx < kA0.getK_eff().size(); ++idx) {
      maxAlphaDiff = std::max(maxAlphaDiff,
                              std::fabs(kAHi.getK_eff()[idx] - kA0.getK_eff()[idx]));
    }
    recordCheck("Physics: CFS alpha_max>0 shifts K_eff", maxAlphaDiff > 1.0, maxAlphaDiff);

    // (d) Kill CPML_NO_CORNER_CROSS: under uniform rigid-body translation u_x = 1,
    //     stiffness force f_K == 0 identically (since sum_a grad N_a == 0), so
    //     1_x^T K_eff 1_x == 1_x^T K_eff^{inertial} 1_x > 0 ONLY if corner cross
    //     terms (g1 o g0) are active!
    double rigidXStiff = 0.0;
    for (int nodeA = 0; nodeA < 4; ++nodeA) {
      for (int nodeB = 0; nodeB < 4; ++nodeB) {
        rigidXStiff += kFull.getK_eff()[(size_t)(nodeA * 2) * nDOF + (nodeB * 2)];
      }
    }
    recordCheck("Physics: corner cross inertial K_{s,mass}>0", rigidXStiff > 1.0, rigidXStiff);
  }
}

static void checkUPPressureShapes(const char* name, CPMLTopology topo,
                                  const std::vector<double>& coords,
                                  double expectedH) {
  const int ndm = shapes::topologyNdm(topo);
  const int nGP = shapes::topologyNumGP(topo);
  const double gotH = shapes::elemSizeH(topo, coords.data());
  const double hErr = std::fabs(gotH - expectedH);

  double worstErr = hErr;
  for (int pOrder = 0; pOrder <= 1; ++pOrder) {
    const int expectedNP = (pOrder == 1) ? shapes::numVertices(topo)
                                         : shapes::topologyNumNodes(topo);
    for (int g = 0; g < nGP; ++g) {
      double Np[shapes::MAX_NODES]{};
      double dNpdx[shapes::MAX_NODES * 3]{};
      int nP_out = -1;
      const bool ok = shapes::evalShapeAndGradP(topo, g, coords.data(),
                                                pOrder, Np, dNpdx, nP_out);
      if (!ok || nP_out != expectedNP) {
        worstErr = 1.0;
        continue;
      }
      double sumNp = 0.0;
      double sumGrad[3] = {0.0, 0.0, 0.0};
      double gradX[3][3]{};
      for (int b = 0; b < nP_out; ++b) {
        sumNp += Np[b];
        for (int j = 0; j < ndm; ++j) {
          sumGrad[j] += dNpdx[b * ndm + j];
          for (int i = 0; i < ndm; ++i) {
            gradX[i][j] += coords[b * ndm + i] * dNpdx[b * ndm + j];
          }
        }
      }
      worstErr = std::max(worstErr, std::fabs(sumNp - 1.0));
      for (int j = 0; j < ndm; ++j) {
        worstErr = std::max(worstErr, std::fabs(sumGrad[j]));
        for (int i = 0; i < ndm; ++i) {
          const double want = (i == j) ? 1.0 : 0.0;
          worstErr = std::max(worstErr, std::fabs(gradX[i][j] - want));
        }
      }
    }
  }
  char label[96];
  std::snprintf(label, sizeof(label), "UP Shapes %s (pOrder=0/1, elemSizeH)", name);
  recordCheck(label, worstErr < 1.0e-12, worstErr);
}

static void checkUPKernelFiniteDiffAndState(const char* label,
                                            CPMLTopology topo,
                                            bool lumped,
                                            const std::vector<double>& coords,
                                            const CPMLGeomConfig& geom,
                                            const CPMLBiotParams& biot) {
  const CPMLMaterialParams mat(2.5e7, 0.25, 2000.0);
  const CPMLGradingParams grading(2.0, 1.0e-6, 3.0, 15.0);

  CPMLElementKernel kernel(topo, lumped, coords.data(), mat, grading, geom, biot);
  if (!kernel.isValid() || !kernel.isUPMode()) {
    recordCheck(label, false, 1.0);
    return;
  }
  const int ndm = kernel.ndm();
  const int nNodes = kernel.numNodes();
  const int nP = kernel.numPressureDOFs();
  const int nDOF = kernel.numDOF();

  // Step 0: advance one step and commit so all Biot recursive convolution memories are active
  const double dt = 0.005;
  std::vector<double> u0(nDOF), v0(nDOF), a0(nDOF);
  for (int i = 0; i < nDOF; ++i) {
    u0[i] = 1.2e-4 * std::sin(0.37 * (i + 1));
    v0[i] = 2.5e-2 * std::cos(0.51 * (i + 1));
    a0[i] = -1.8e0 * std::sin(0.29 * (i + 1));
  }
  kernel.setTrial(u0.data(), v0.data(), a0.data(), dt);
  kernel.commit();

  std::vector<double> packedAfterCommit(kernel.numStateDoubles(), 0.0);
  kernel.packState(packedAfterCommit.data());

  // Step 1: evaluate at a second non-trivial trial state
  std::vector<double> u1(nDOF), v1(nDOF), a1(nDOF);
  for (int i = 0; i < nDOF; ++i) {
    u1[i] = u0[i] + 0.8e-4 * std::cos(0.43 * (i + 1));
    v1[i] = v0[i] - 1.4e-2 * std::sin(0.61 * (i + 1));
    a1[i] = a0[i] + 2.1e0  * std::cos(0.17 * (i + 1));
  }
  kernel.setTrial(u1.data(), v1.data(), a1.data(), dt);

  const double c1 = 1.0;
  const double c2 = 2.0 / dt;
  const double c3 = 4.0 / (dt * dt);
  const std::vector<double> R_base = kernel.getR_int();
  const std::vector<double> K_eff  = kernel.getK_eff();
  const std::vector<double> C_eff  = kernel.getC_eff();
  const std::vector<double> M_eff  = kernel.getM_eff();
  const std::vector<double> KT_dyn = kernel.getK_T_dyn(c1, c2, c3);

  // Verify skew-duality K_eff[uRow, pCol] == -C_eff[pCol, uRow] and zero blocks
  double maxScale = 1.0;
  for (double val : K_eff) maxScale = std::max(maxScale, std::fabs(val));
  for (double val : C_eff) maxScale = std::max(maxScale, std::fabs(val));
  double skewAndBlockErr = 0.0;
  for (int nodeA = 0; nodeA < nNodes; ++nodeA) {
    for (int d = 0; d < ndm; ++d) {
      const int uRow = kernel.uOffset(nodeA) + d;
      for (int b = 0; b < nP; ++b) {
        const int pCol = kernel.pOffset(b);
        const double k_up = K_eff[(size_t)uRow * nDOF + pCol];
        const double c_pu = C_eff[(size_t)pCol * nDOF + uRow];
        const double k_pu = K_eff[(size_t)pCol * nDOF + uRow];
        const double c_up = C_eff[(size_t)uRow * nDOF + pCol];
        const double m_up = M_eff[(size_t)uRow * nDOF + pCol];
        const double m_pu = M_eff[(size_t)pCol * nDOF + uRow];
        skewAndBlockErr = std::max(skewAndBlockErr, std::fabs(k_up + c_pu) / maxScale);
        skewAndBlockErr = std::max(skewAndBlockErr, std::fabs(k_pu) / maxScale);
        skewAndBlockErr = std::max(skewAndBlockErr, std::fabs(c_up) / maxScale);
        skewAndBlockErr = std::max(skewAndBlockErr, std::fabs(m_up) / maxScale);
        skewAndBlockErr = std::max(skewAndBlockErr, std::fabs(m_pu) / maxScale);
      }
    }
  }
  for (int b = 0; b < nP; ++b) {
    const int pRow = kernel.pOffset(b);
    for (int c = 0; c < nP; ++c) {
      const int pCol = kernel.pOffset(c);
      skewAndBlockErr = std::max(skewAndBlockErr, std::fabs(M_eff[(size_t)pRow * nDOF + pCol]));
      skewAndBlockErr = std::max(skewAndBlockErr,
          std::fabs(K_eff[(size_t)pRow * nDOF + pCol] - K_eff[(size_t)pCol * nDOF + pRow]) / maxScale);
      skewAndBlockErr = std::max(skewAndBlockErr,
          std::fabs(C_eff[(size_t)pRow * nDOF + pCol] - C_eff[(size_t)pCol * nDOF + pRow]) / maxScale);
    }
  }
  char skewMsg[96];
  std::snprintf(skewMsg, sizeof(skewMsg), "UP Skew-duality & blocks %s", label);
  recordCheck(skewMsg, skewAndBlockErr < 1.0e-13, skewAndBlockErr);

  // Verify partial FD tangents dR/dU = K_eff, dR/dV = C_eff, dR/dA = M_eff, and Newmark K_T_dyn
  // using per-block scales (uu, up, pu, pp) so O(1e-8) Darcy/storage blocks are verified!
  const double eps = 1.0e-7;
  double maxFdErr = 0.0;
  double scK[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  double scC[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  double scM[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  double scKT[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  for (int i = 0; i < nDOF; ++i) {
    const int bi = isPressureDOF(kernel, i) ? 1 : 0;
    for (int j = 0; j < nDOF; ++j) {
      const int bj = isPressureDOF(kernel, j) ? 1 : 0;
      const size_t idx = (size_t)i * nDOF + j;
      scK[bi][bj]  = std::max(scK[bi][bj],  std::fabs(K_eff[idx]));
      scC[bi][bj]  = std::max(scC[bi][bj],  std::fabs(C_eff[idx]));
      scM[bi][bj]  = std::max(scM[bi][bj],  std::fabs(M_eff[idx]));
      scKT[bi][bj] = std::max(scKT[bi][bj], std::fabs(KT_dyn[idx]));
    }
  }
  for (int bi = 0; bi < 2; ++bi) {
    for (int bj = 0; bj < 2; ++bj) {
      if (scK[bi][bj]  <= 0.0) scK[bi][bj]  = 1.0;
      if (scC[bi][bj]  <= 0.0) scC[bi][bj]  = 1.0;
      if (scM[bi][bj]  <= 0.0) scM[bi][bj]  = 1.0;
      if (scKT[bi][bj] <= 0.0) scKT[bi][bj] = 1.0;
    }
  }

  for (int j = 0; j < nDOF; ++j) {
    const int bj = isPressureDOF(kernel, j) ? 1 : 0;
    const double eps_j = bj ? (eps * mat.E) : eps;
    // 1. dR_disp/dU_j == K_eff[:, j] (covers both dR/du and dR/dp across all 4 blocks)
    {
      std::vector<double> up = u1, um = u1;
      up[j] += eps_j; um[j] -= eps_j;
      kernel.setTrial(up.data(), v1.data(), a1.data(), dt);
      const std::vector<double> Rp = kernel.getR_disp();
      kernel.setTrial(um.data(), v1.data(), a1.data(), dt);
      const std::vector<double> Rm = kernel.getR_disp();
      for (int i = 0; i < nDOF; ++i) {
        const int bi = isPressureDOF(kernel, i) ? 1 : 0;
        const double fd = (Rp[i] - Rm[i]) / (2.0 * eps_j);
        maxFdErr = std::max(maxFdErr, std::fabs(fd - K_eff[(size_t)i * nDOF + j]) / scK[bi][bj]);
      }
    }
    // 2. dR_vel/dV_j == C_eff[:, j] (covers both dR/dv_u and dR/dv_p across all 4 blocks)
    {
      std::vector<double> vp = v1, vm = v1;
      vp[j] += eps_j; vm[j] -= eps_j;
      kernel.setTrial(u1.data(), vp.data(), a1.data(), dt);
      const std::vector<double> Rp = kernel.getR_vel();
      kernel.setTrial(u1.data(), vm.data(), a1.data(), dt);
      const std::vector<double> Rm = kernel.getR_vel();
      for (int i = 0; i < nDOF; ++i) {
        const int bi = isPressureDOF(kernel, i) ? 1 : 0;
        const double fd = (Rp[i] - Rm[i]) / (2.0 * eps_j);
        maxFdErr = std::max(maxFdErr, std::fabs(fd - C_eff[(size_t)i * nDOF + j]) / scC[bi][bj]);
      }
    }
    // 3. dR_int/dA_j == M_eff[:, j]
    {
      const double eps_a = c3 * eps;
      std::vector<double> ap = a1, am = a1;
      ap[j] += eps_a; am[j] -= eps_a;
      kernel.setTrial(u1.data(), v1.data(), ap.data(), dt);
      const std::vector<double> Rp = kernel.getR_int();
      kernel.setTrial(u1.data(), v1.data(), am.data(), dt);
      const std::vector<double> Rm = kernel.getR_int();
      for (int i = 0; i < nDOF; ++i) {
        const int bi = isPressureDOF(kernel, i) ? 1 : 0;
        const double fd = (Rp[i] - Rm[i]) / (2.0 * eps_a);
        maxFdErr = std::max(maxFdErr, std::fabs(fd - M_eff[(size_t)i * nDOF + j]) / scM[bi][bj]);
      }
    }
    // 4. Combined Newmark K_T_dyn = c1*dR_int/dU + c2*dR_int/dV + c3*dR_int/dA
    {
      std::vector<double> up = u1, vp = v1, ap = a1;
      std::vector<double> um = u1, vm = v1, am = a1;
      up[j] += c1 * eps_j; vp[j] += c2 * eps_j; ap[j] += c3 * eps_j;
      um[j] -= c1 * eps_j; vm[j] -= c2 * eps_j; am[j] -= c3 * eps_j;
      kernel.setTrial(up.data(), vp.data(), ap.data(), dt);
      const std::vector<double> Rp = kernel.getR_int();
      kernel.setTrial(um.data(), vm.data(), am.data(), dt);
      const std::vector<double> Rm = kernel.getR_int();
      for (int i = 0; i < nDOF; ++i) {
        const int bi = isPressureDOF(kernel, i) ? 1 : 0;
        const double fd = (Rp[i] - Rm[i]) / (2.0 * eps_j);
        maxFdErr = std::max(maxFdErr, std::fabs(fd - KT_dyn[(size_t)i * nDOF + j]) / scKT[bi][bj]);
      }
    }
  }
  char fdMsg[96];
  std::snprintf(fdMsg, sizeof(fdMsg), "UP 4-Block & Newmark FD %s", label);
  recordCheck(fdMsg, maxFdErr < 1.0e-6, maxFdErr);

  // Verify revertToLastCommit, packState/unpackState, and revertToStart in -up mode
  kernel.revertToLastCommit();
  std::vector<double> packedAfterRevert(kernel.numStateDoubles(), 0.0);
  kernel.packState(packedAfterRevert.data());
  double revertPackDiff = 0.0;
  for (size_t i = 0; i < packedAfterCommit.size(); ++i) {
    revertPackDiff = std::max(revertPackDiff,
                              std::fabs(packedAfterRevert[i] - packedAfterCommit[i]));
  }
  kernel.setTrial(u1.data(), v1.data(), a1.data(), dt);
  const std::vector<double> R_re = kernel.getR_int();
  double revErr = 0.0;
  for (int i = 0; i < nDOF; ++i) {
    revErr = std::max(revErr, std::fabs(R_re[i] - R_base[i]));
  }

  CPMLElementKernel clone(topo, lumped, coords.data(), mat, grading, geom, biot);
  clone.unpackState(packedAfterCommit.data());
  clone.setTrial(u1.data(), v1.data(), a1.data(), dt);
  double packErr = 0.0;
  for (int i = 0; i < nDOF; ++i) {
    packErr = std::max(packErr, std::fabs(clone.getR_int()[i] - R_base[i]));
  }

  kernel.revertToStart();
  kernel.setTrial(u0.data(), v0.data(), a0.data(), dt);
  CPMLElementKernel fresh(topo, lumped, coords.data(), mat, grading, geom, biot);
  fresh.setTrial(u0.data(), v0.data(), a0.data(), dt);
  double startErr = 0.0;
  for (int i = 0; i < nDOF; ++i) {
    startErr = std::max(startErr, std::fabs(kernel.getR_int()[i] - fresh.getR_int()[i]));
  }

  char stateMsg[96];
  std::snprintf(stateMsg, sizeof(stateMsg), "UP State/Revert/Pack %s", label);
  const double worstStateErr = std::max({revertPackDiff, revErr, packErr, startErr});
  recordCheck(stateMsg, worstStateErr == 0.0, worstStateErr);
}

static void checkUPPhysicsAndStaticLimits(const std::vector<double>& xyQ4) {
  const CPMLMaterialParams mat(2.5e7, 0.25, 2000.0);
  const CPMLGradingParams grading(2.0, 1.0e-6, 3.0, 15.0);
  const CPMLGeomConfig cornerBox = CPMLGeomConfig::makeBox2D(-1.0, 0.05, -1.0, 0.08, 1.5, 1.2);
  const CPMLGeomConfig interiorBox = CPMLGeomConfig::makeBox2D(-10.0, 10.0, -10.0, 10.0, 1.5, 1.2);

  CPMLBiotParams biot;
  biot.enabled = true;
  biot.Kf = 2.2e9;
  biot.poro = 0.35;
  biot.rhoF = 1000.0;
  biot.biotAlpha = 0.85;
  biot.Ks = 3.6e10;
  biot.pOrder = 0;
  biot.formulation = 1; // bbar
  biot.stabMode = 1;    // auto
  biot.stabValue = 0.25;
  biot.kbar[0] = 1.0e-7;
  biot.kbar[1] = 2.0e-7;

  // 1. Undrained Vp > drained Vp
  const double vpU = biot.undrainedVp(mat);
  recordCheck("UP Physics: undrained Vp > drained Vp", vpU > mat.Vp, vpU - mat.Vp);

  // 2. Interior GP limit in -up mode: K_eff == K_0, C_eff == C_0, M_eff == M_0, and forces
  {
    CPMLElementKernel kInt(CPMLTopology::Q4, false, xyQ4.data(), mat, grading, interiorBox, biot);
    const int nDOF = kInt.numDOF();
    std::vector<double> u(nDOF), v(nDOF), a(nDOF);
    for (int i = 0; i < nDOF; ++i) {
      u[i] = 1.0e-3 * (i + 1);
      v[i] = 0.1 * (i + 1);
      a[i] = 0.5 * (i + 1);
    }
    kInt.setTrial(u.data(), v.data(), a.data(), 0.005);
    double intErr = 0.0;
    for (size_t idx = 0; idx < kInt.getK0().size(); ++idx) {
      intErr = std::max(intErr, std::fabs(kInt.getK_eff()[idx] - kInt.getK0()[idx]));
      intErr = std::max(intErr, std::fabs(kInt.getC_eff()[idx] - kInt.getC0()[idx]));
      intErr = std::max(intErr, std::fabs(kInt.getM_eff()[idx] - kInt.getM0()[idx]));
    }
    for (int i = 0; i < nDOF; ++i) {
      double wantDisp = 0.0, wantVel = 0.0, wantAcc = 0.0;
      for (int j = 0; j < nDOF; ++j) {
        const size_t idx = (size_t)i * nDOF + j;
        wantDisp += kInt.getK0()[idx] * u[j];
        wantVel  += kInt.getC0()[idx] * v[j];
        wantAcc  += kInt.getM0()[idx] * a[j];
      }
      const double sc = std::max(1.0, std::fabs(wantDisp) + std::fabs(wantVel) + std::fabs(wantAcc));
      intErr = std::max(intErr, std::fabs(kInt.getR_disp()[i] - wantDisp) / sc);
      intErr = std::max(intErr, std::fabs(kInt.getR_vel()[i] - wantVel) / sc);
      intErr = std::max(intErr, std::fabs(kInt.getR_int()[i] - (wantDisp + wantVel + wantAcc)) / sc);
    }
    recordCheck("UP Interior GP limit (xi=0 -> K0, C0, M0 & forces)", intErr < 1.0e-12, intErr);
  }

  // 3. Static limit (dt <= 0.0 and setStaticMode(true)) in -up mode
  {
    CPMLElementKernel kStat(CPMLTopology::Q4, false, xyQ4.data(), mat, grading, cornerBox, biot);
    const int nDOF = kStat.numDOF();
    std::vector<double> u(nDOF), v(nDOF), a(nDOF);
    for (int i = 0; i < nDOF; ++i) {
      u[i] = 1.0e-3 * (i + 1);
      v[i] = 0.2 * (i + 1);
      a[i] = 1.5 * (i + 1);
    }
    kStat.setTrial(u.data(), v.data(), a.data(), 0.005);
    kStat.commit();
    std::vector<double> stateBeforeStatic(kStat.numStateDoubles(), 0.0);
    kStat.packState(stateBeforeStatic.data());

    kStat.setTrial(u.data(), v.data(), a.data(), 0.0);
    double statErr = 0.0;
    for (size_t idx = 0; idx < kStat.getK0().size(); ++idx) {
      statErr = std::max(statErr, std::fabs(kStat.getK_eff()[idx] - kStat.getK0()[idx]));
      statErr = std::max(statErr, std::fabs(kStat.getC_eff()[idx] - kStat.getC0()[idx]));
    }
    for (int i = 0; i < nDOF; ++i) {
      double wantDisp = 0.0, wantVel = 0.0, wantAcc = 0.0;
      for (int j = 0; j < nDOF; ++j) {
        const size_t idx = (size_t)i * nDOF + j;
        wantDisp += kStat.getK0()[idx] * u[j];
        wantVel  += kStat.getC0()[idx] * v[j];
        wantAcc  += kStat.getM0()[idx] * a[j];
      }
      statErr = std::max(statErr, std::fabs(kStat.getR_disp()[i] - wantDisp));
      statErr = std::max(statErr, std::fabs(kStat.getR_vel()[i] - wantVel));
      statErr = std::max(statErr, std::fabs(kStat.getR_int()[i] - (wantDisp + wantVel + wantAcc)));
    }
    kStat.setStaticMode(true);
    kStat.setTrial(u.data(), v.data(), a.data(), 0.005);
    for (size_t idx = 0; idx < kStat.getK0().size(); ++idx) {
      statErr = std::max(statErr, std::fabs(kStat.getK_eff()[idx] - kStat.getK0()[idx]));
      statErr = std::max(statErr, std::fabs(kStat.getC_eff()[idx] - kStat.getC0()[idx]));
    }
    kStat.setStaticMode(false);
    kStat.revertToLastCommit();
    std::vector<double> stateAfterStatic(kStat.numStateDoubles(), 0.0);
    kStat.packState(stateAfterStatic.data());
    for (size_t idx = 0; idx < stateBeforeStatic.size(); ++idx) {
      statErr = std::max(statErr, std::fabs(stateAfterStatic[idx] - stateBeforeStatic[idx]));
    }
    recordCheck("UP Static limit (dt<=0 & staticMode -> K0, C0, M0)", statErr == 0.0, statErr);
  }

  // 4. CPMLBiotParams::isValid negative validation & degenerate (1/Qbar=0, kbar=0) checks
  {
    auto badPoroLow = biot; badPoroLow.poro = -0.01;
    auto badPoroOne = biot; badPoroOne.poro = 1.0;
    auto badKf      = biot; badKf.Kf = 0.0;
    auto badRhoF    = biot; badRhoF.rhoF = -1.0;
    auto badAlphaLo = biot; badAlphaLo.biotAlpha = 0.20; // < poro=0.35
    auto badAlphaHi = biot; badAlphaHi.biotAlpha = 1.05;
    auto badKs      = biot; badKs.Ks = -1.0e9;
    auto badKbar    = biot; badKbar.kbar[0] = -1.0e-8;
    auto badPOrder  = biot; badPOrder.pOrder = 2;
    auto badForm    = biot; badForm.formulation = 2;
    auto badStabM   = biot; badStabM.stabMode = 3;
    auto badStabV   = biot; badStabV.stabMode = 0; badStabV.stabValue = -0.1;

    const bool validChecks =
        biot.isValid(2, CPMLTopology::Q4) &&
        !badPoroLow.isValid(2, CPMLTopology::Q4) &&
        !badPoroOne.isValid(2, CPMLTopology::Q4) &&
        !badKf.isValid(2, CPMLTopology::Q4) &&
        !badRhoF.isValid(2, CPMLTopology::Q4) &&
        !badAlphaLo.isValid(2, CPMLTopology::Q4) &&
        !badAlphaHi.isValid(2, CPMLTopology::Q4) &&
        !badKs.isValid(2, CPMLTopology::Q4) &&
        !badKbar.isValid(2, CPMLTopology::Q4) &&
        !badPOrder.isValid(2, CPMLTopology::Q4) &&
        !badForm.isValid(2, CPMLTopology::Q4) &&
        !badStabM.isValid(2, CPMLTopology::Q4) &&
        !badStabV.isValid(2, CPMLTopology::Q4);

    // Degenerate 1/Qbar = 0 and kbar = 0 (impermeable, incompressible limit)
    auto degen = biot;
    degen.poro = 0.0;
    degen.Ks = 0.0;
    degen.kbar[0] = 0.0;
    degen.kbar[1] = 0.0;
    CPMLElementKernel kDeg(CPMLTopology::Q4, false, xyQ4.data(), mat, grading, cornerBox, degen);
    std::vector<double> u(kDeg.numDOF(), 1.0e-4), v(kDeg.numDOF(), 0.01), a(kDeg.numDOF(), 0.1);
    kDeg.setTrial(u.data(), v.data(), a.data(), 0.002);
    bool finiteDeg = kDeg.isValid();
    for (double val : kDeg.getR_int()) {
      if (!std::isfinite(val)) finiteDeg = false;
    }
    recordCheck("UP BiotParams::isValid & degenerate 1/Qbar=0, kbar=0", validChecks && finiteDeg);
  }
}

static void runSelfChecks() {
  std::printf("=== RUNNING LADRUNOCPML C++ SELF-CHECKS ===\n");

  // Canonical positively-oriented element coordinates for all 9 topologies
  // 1. T3: (0,0), (1,0), (0,1) -> area = 0.5
  const std::vector<double> xyT3 = {
    0.0, 0.0,  1.0, 0.0,  0.0, 1.0
  };
  // 2. Q4: [0,1]x[0,1] -> area = 1.0
  const std::vector<double> xyQ4 = {
    0.0, 0.0,  1.0, 0.0,  1.0, 1.0,  0.0, 1.0
  };
  // 3. T6: Node 0=(1,0), 1=(0,1), 2=(0,0), 3=(0.5,0.5), 4=(0,0.5), 5=(0.5,0) -> area = 0.5
  const std::vector<double> xyT6 = {
    1.0, 0.0,  0.0, 1.0,  0.0, 0.0,
    0.5, 0.5,  0.0, 0.5,  0.5, 0.0
  };
  // 4. BT6: Node 0=(0,0), 1=(1,0), 2=(0,1), 3=(0.5,0), 4=(0.5,0.5), 5=(0,0.5) -> area = 0.5
  const std::vector<double> xyBT6 = {
    0.0, 0.0,  1.0, 0.0,  0.0, 1.0,
    0.5, 0.0,  0.5, 0.5,  0.0, 0.5
  };
  // 5. TET4: Node 0=(0,0,0), 1=(1,0,0), 2=(0,1,0), 3=(0,0,1) -> vol = 1/6
  const std::vector<double> xyzTET4 = {
    0.0, 0.0, 0.0,  1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0
  };
  // 6. H8: [0,1]^3 -> vol = 1.0
  const std::vector<double> xyzH8 = {
    0.0, 0.0, 0.0,  1.0, 0.0, 0.0,  1.0, 1.0, 0.0,  0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,  1.0, 0.0, 1.0,  1.0, 1.0, 1.0,  0.0, 1.0, 1.0
  };
  // 7. TET10 & 8. BTET10 (detJ > 0 ordering): Node 0=(1,0,0), 1=(0,1,0), 2=(0,0,1), 3=(0,0,0)
  const std::vector<double> xyzTET10 = {
    1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0,  0.0, 0.0, 0.0,
    0.5, 0.5, 0.0,  0.0, 0.5, 0.5,  0.5, 0.0, 0.5,
    0.5, 0.0, 0.0,  0.0, 0.0, 0.5,  0.0, 0.5, 0.0
  };
  // Standard right-handed OpenSees TenNodeTetrahedron / BezierTet10 ordering (detJ < 0):
  // Node 0=(0,0,0), 1=(1,0,0), 2=(0,1,0), 3=(0,0,1) + straight-edge midpoints
  const std::vector<double> xyzTET10_RH = {
    0.0, 0.0, 0.0,  1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0,
    0.5, 0.0, 0.0,  0.5, 0.5, 0.0,  0.0, 0.5, 0.0,
    0.0, 0.0, 0.5,  0.0, 0.5, 0.5,  0.5, 0.0, 0.5
  };
  // 9. H20: [0,1]^3 -> vol = 1.0
  const std::vector<double> xyzH20 = {
    0.0, 0.0, 0.0,  1.0, 0.0, 0.0,  1.0, 1.0, 0.0,  0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,  1.0, 0.0, 1.0,  1.0, 1.0, 1.0,  0.0, 1.0, 1.0,
    0.5, 0.0, 0.0,  1.0, 0.5, 0.0,  0.5, 1.0, 0.0,  0.0, 0.5, 0.0,
    0.5, 0.0, 1.0,  1.0, 0.5, 1.0,  0.5, 1.0, 1.0,  0.0, 0.5, 1.0,
    0.0, 0.0, 0.5,  1.0, 0.0, 0.5,  1.0, 1.0, 0.5,  0.0, 1.0, 0.5
  };

  checkTopologyShapes("T3",        CPMLTopology::T3,     xyT3,        0.5);
  checkTopologyShapes("Q4",        CPMLTopology::Q4,     xyQ4,        1.0);
  checkTopologyShapes("T6",        CPMLTopology::T6,     xyT6,        0.5);
  checkTopologyShapes("BT6",       CPMLTopology::BT6,    xyBT6,       0.5);
  checkTopologyShapes("TET4",      CPMLTopology::TET4,   xyzTET4,     1.0 / 6.0);
  checkTopologyShapes("H8",        CPMLTopology::H8,     xyzH8,       1.0);
  checkTopologyShapes("TET10",     CPMLTopology::TET10,  xyzTET10,    1.0 / 6.0);
  checkTopologyShapes("BTET10",    CPMLTopology::BTET10, xyzTET10,    1.0 / 6.0);
  checkTopologyShapes("TET10-RH",  CPMLTopology::TET10,  xyzTET10_RH, 1.0 / 6.0, true);
  checkTopologyShapes("BTET10-RH", CPMLTopology::BTET10, xyzTET10_RH, 1.0 / 6.0, true);
  checkTopologyShapes("H20",       CPMLTopology::H20,    xyzH20,      1.0);

  const double sqrt2 = std::sqrt(2.0);
  checkUPPressureShapes("T3",     CPMLTopology::T3,     xyT3,        sqrt2);
  checkUPPressureShapes("Q4",     CPMLTopology::Q4,     xyQ4,        1.0);
  checkUPPressureShapes("T6",     CPMLTopology::T6,     xyT6,        sqrt2);
  checkUPPressureShapes("BT6",    CPMLTopology::BT6,    xyBT6,       sqrt2);
  checkUPPressureShapes("TET4",   CPMLTopology::TET4,   xyzTET4,     sqrt2);
  checkUPPressureShapes("H8",     CPMLTopology::H8,     xyzH8,       1.0);
  checkUPPressureShapes("TET10",  CPMLTopology::TET10,  xyzTET10,    sqrt2);
  checkUPPressureShapes("BTET10", CPMLTopology::BTET10, xyzTET10_RH, sqrt2);
  checkUPPressureShapes("H20",    CPMLTopology::H20,    xyzH20,      1.0);

  // Corner-activated 2D/3D Box and oblique Plane geometries
  const CPMLGeomConfig box2D   = CPMLGeomConfig::makeBox2D(-1.0, 0.05, -1.0, 0.08, 1.5, 1.2);
  const CPMLGeomConfig plane2D = CPMLGeomConfig::makePlane2D(-0.1, -0.1, 0.6, 0.8, 1.5);
  const CPMLGeomConfig box3D   = CPMLGeomConfig::makeBox3D(-1.0, 0.05, -1.0, 0.08, -1.0, 0.06, 1.5, 1.2, 1.4);
  const CPMLGeomConfig plane3D = CPMLGeomConfig::makePlane3D(-0.1, -0.1, -0.1, 0.3, 0.4, 0.8660254037844386, 1.5);

  checkKernelFiniteDiffAndState("T3-box-cons",          CPMLTopology::T3,     false, xyT3,        box2D);
  checkKernelFiniteDiffAndState("Q4-plane-lump",        CPMLTopology::Q4,     true,  xyQ4,        plane2D);
  checkKernelFiniteDiffAndState("T6-box-lump",          CPMLTopology::T6,     true,  xyT6,        box2D);
  checkKernelFiniteDiffAndState("BT6-plane-cons",       CPMLTopology::BT6,    false, xyBT6,       plane2D);
  checkKernelFiniteDiffAndState("TET4-box-cons",        CPMLTopology::TET4,   false, xyzTET4,     box3D);
  checkKernelFiniteDiffAndState("H8-plane-lump",        CPMLTopology::H8,     true,  xyzH8,       plane3D);
  checkKernelFiniteDiffAndState("TET10-box-lump",       CPMLTopology::TET10,  true,  xyzTET10,    box3D);
  checkKernelFiniteDiffAndState("TET10-RH-box-cons",    CPMLTopology::TET10,  false, xyzTET10_RH, box3D);
  checkKernelFiniteDiffAndState("BTET10-plane-cons",    CPMLTopology::BTET10, false, xyzTET10,    plane3D);
  checkKernelFiniteDiffAndState("BTET10-RH-plane-lump", CPMLTopology::BTET10, true,  xyzTET10_RH, plane3D);
  checkKernelFiniteDiffAndState("H20-box-lump",         CPMLTopology::H20,    true,  xyzH20,      box3D);

  // Biot u-p kernel self-checks across equal-order (pOrder=0) and Taylor-Hood (pOrder=1),
  // std (formulation=0) and bbar (formulation=1), and stabilization modes
  CPMLBiotParams biotEqStd;
  biotEqStd.enabled = true;
  biotEqStd.Kf = 2.2e9;
  biotEqStd.poro = 0.35;
  biotEqStd.rhoF = 1000.0;
  biotEqStd.biotAlpha = 0.85;
  biotEqStd.Ks = 3.6e10;
  biotEqStd.pOrder = 0;
  biotEqStd.formulation = 0;
  biotEqStd.stabMode = 1;
  biotEqStd.stabValue = 0.25;
  biotEqStd.kbar[0] = 1.2e-7;
  biotEqStd.kbar[1] = 2.5e-7;
  biotEqStd.kbar[2] = 0.8e-7;

  CPMLBiotParams biotEqBbar = biotEqStd;
  biotEqBbar.formulation = 1;

  CPMLBiotParams biotTHStd = biotEqStd;
  biotTHStd.pOrder = 1;
  biotTHStd.formulation = 0;
  biotTHStd.stabMode = 0;

  CPMLBiotParams biotTHBbar = biotTHStd;
  biotTHBbar.formulation = 1;

  checkUPKernelFiniteDiffAndState("T3-UP-eq-std-box",     CPMLTopology::T3,     false, xyT3,        box2D,   biotEqStd);
  checkUPKernelFiniteDiffAndState("Q4-UP-eq-bbar-plane",  CPMLTopology::Q4,     true,  xyQ4,        plane2D, biotEqBbar);
  checkUPKernelFiniteDiffAndState("T6-UP-TH-std-box",     CPMLTopology::T6,     false, xyT6,        box2D,   biotTHStd);
  checkUPKernelFiniteDiffAndState("BT6-UP-TH-bbar-plane", CPMLTopology::BT6,    true,  xyBT6,       plane2D, biotTHBbar);
  checkUPKernelFiniteDiffAndState("TET4-UP-eq-std-box",   CPMLTopology::TET4,   false, xyzTET4,     box3D,   biotEqStd);
  checkUPKernelFiniteDiffAndState("H8-UP-eq-bbar-plane",  CPMLTopology::H8,     true,  xyzH8,       plane3D, biotEqBbar);
  checkUPKernelFiniteDiffAndState("TET10-UP-TH-std-box",  CPMLTopology::TET10,  false, xyzTET10_RH, box3D,   biotTHStd);
  checkUPKernelFiniteDiffAndState("BTET10-UP-TH-bbar-pl", CPMLTopology::BTET10, true,  xyzTET10_RH, plane3D, biotTHBbar);
  checkUPKernelFiniteDiffAndState("H20-UP-TH-bbar-box",   CPMLTopology::H20,    false, xyzH20,      box3D,   biotTHBbar);

  checkPhysicsLimitsAndMutations(xyQ4);
  checkUPPhysicsAndStaticLimits(xyQ4);
}

int main(int argc, char** argv) {
  bool selfCheckOnly = false;
  std::vector<std::string> casePaths;

  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--self-check") == 0) {
      selfCheckOnly = true;
    } else {
      casePaths.emplace_back(argv[i]);
    }
  }
  if (!selfCheckOnly && casePaths.empty()) {
    casePaths.emplace_back("tests/ladruno_cpml_cases.txt");
  }

  runSelfChecks();

  if (!selfCheckOnly) {
    for (const auto& path : casePaths) {
      int nCases = runCasesFile(path);
      if (nCases <= 0) {
        g_failures++;
      }
    }
  }

  if (g_failures > 0) {
    std::printf("\nRESULT: FAIL (%d failure(s))\n", g_failures);
    return 1;
  }
  std::printf("\nRESULT: PASS\n");
  return 0;
}
