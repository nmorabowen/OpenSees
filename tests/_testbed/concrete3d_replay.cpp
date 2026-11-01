// concrete3d_replay.cpp -- record-replay testbed for LadrunoConcrete3DKernel.h (apeGmsh #2621, increment F1).
//
// The fixture (concrete3d_replay_records.txt) is a subset of the LADRUNO_C3D_DIAG record dump written by the shipped
// opensees.pyd on the #2602 corbel (fork probe/2621-recdump, a counters-only instrumentation of d14fb6958). A record
// carries EVERYTHING the kernel reads (every Params field, the full committed State, the trial strain) and what the
// binary returned. This testbed puts each record back through the kernel OF THIS TREE and requires the recorded answer
// BIT FOR BIT; nothing is compared with a tolerance and no field or record class is exempt.
//   C, F  one integrate(): K::returnMap exactly as LadrunoConcrete3D::integrate calls it. Compared: status, piece count
//         (subInfo), kappa, omega_t, omega_c, nominal and effective stress, the full 6x6 tangent, and eig3sym of the
//         effective stress (W line: the libm canary, eig3sym is atan2/cos/acos-heavy).
//   R     one refused plane-stress setTrialStrain: the sigma_22 Newton of LadrunoConcrete3D::setTrialStrain (rule B,
//         copied below) from the recorded eps_22 start. Compared: the refusal reason, the pass count and every pass's
//         eps_22, sigma_22, D_22, status and piece count, plus the last pass's outcome (O), full 6x6 tangent and eig3sym.
// Tangent lines: D = written by the pyd (it dumped the 6x6 only on C records); M = the 6x6 the pyd computed but did not
// dump (F records, the last pass of R records), recorded by this testbed under the same libm on the same kernel blob,
// after every pyd-written value of the record replayed bit-exactly. A C record must carry D, an R or F record M.
// Bit-exactness holds only on the libm the records were written with: Intel's libmmt for the Windows pyd (its Fortran
// objects pull it in by /DEFAULTLIB), glibc for the Linux build. The build defines C3D_REPLAY_LIBM to the libm it
// links (for glibc linux-glibc-<the version ldd --version reports>, checked against the loaded libc at start-up), and
// `check` refuses a fixture recorded under another one, another glibc version included (--cross compares anyway).
// A KERNEL CHANGE THAT MOVES RESULTS FAILS `check` BY DESIGN. Such a PR re-records the fixture on each libm (`record`)
// and commits it; the fixture's diff is the exact list of values the change moved, and the reviewers judge that diff.
//
//   concrete3d_replay check  <fixture> [--cross]                 exit 0 iff every record replays bit-exactly
//   concrete3d_replay record <fixture> <out> --kernel-blob <sha> re-record every output on THIS build (inputs, order
//                                                                and selection kept), stamped with this libm and the
//                                                                kernel header's git blob
// Builds: concrete3d_replay_ci.sh (g++, glibc) and concrete3d_replay_libmmt.bat (MSVC + Intel libmmt, as the pyd).
// The sigma_22 loop below is a COPY of the wrapper's: the fixture's R records prove the copy matches the wrapper the
// records were written by; a later change to LadrunoConcrete3D::setTrialStrain must be mirrored here.
#include "SRC/material/nD/LadrunoConcrete3DKernel.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#ifndef C3D_REPLAY_LIBM
#error "define C3D_REPLAY_LIBM to the libm this build links (win-libmmt, win-ucrt, linux-glibc-<ldd --version>)"
#endif
#ifdef __GLIBC__
#include <gnu/libc-version.h>
#endif

namespace K = Ladruno::Concrete3D;

[[noreturn]] static void die(const std::string& m) { std::fprintf(stderr, "REPLAY ERROR: %s\n", m.c_str()); std::exit(2); }

struct Pass { int p, status, sub; double eps22, s22, d22; };
struct Rec {
    std::string kLine, gLine; std::map<std::string, std::string> k;   // K line (verbatim + parsed), G line (selection)
    char kind = 0; int code = 0;
    std::map<std::string, std::string> pkv; std::string pLine;
    std::vector<std::string> stateLines;                              // E S X H I C T, verbatim
    double eps[6], sig[6], sigEff[6], h[11], im[11], c[2], strain[6];
    bool hasO = false, hasD = false, hasW = false; char dKind = 0;   // 'D' (pyd-written) or 'M' (testbed-recorded)
    double o[17], D[36], w[12];
    std::vector<Pass> z;
};

static std::vector<double> nums(const std::string& line, size_t want) {
    std::istringstream is(line.substr(1)); std::vector<double> v; std::string t;
    while (is >> t) { char* end = 0; const double x = std::strtod(t.c_str(), &end); if (*end) die("not a number: " + t + " in: " + line.substr(0, 60)); v.push_back(x); }
    if (want && v.size() != want) die("line has " + std::to_string(v.size()) + " values, expected " + std::to_string(want) + ": " + line.substr(0, 60));
    return v;
}
static void kv(const std::string& line, std::map<std::string, std::string>& m) {
    std::istringstream is(line.substr(1)); std::string t;
    while (is >> t) { const size_t p = t.find('='); if (p == std::string::npos) die("not key=value: " + t); m[t.substr(0, p)] = t.substr(p + 1); }
}

// An integer field (status, piece count, pass index) parsed as a double: a non-integral or non-finite value is an error,
// never truncated (a truncating cast would let a corrupted 0.9 replay as 0).
static int asInt(double v, const char* what, const std::string& line) {
    if (!std::isfinite(v) || v != std::trunc(v)) die(std::string(what) + " is not an integer: " + line.substr(0, 60));
    return (int)v;
}

struct Fixture { std::map<std::string, std::string> meta; std::string qLine; std::vector<Rec> recs; std::vector<std::string> header; };

static Fixture load(const char* path) {
    std::ifstream in(path); if (!in) die(std::string("cannot open ") + path);
    Fixture fx; std::string line; Rec* r = 0;
    auto copy = [](const std::vector<double>& v, double* d) { for (size_t i = 0; i < v.size(); ++i) d[i] = v[i]; };
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '#') { if (line.rfind("#@ ", 0) == 0) kv(line.substr(1), fx.meta); if (!r) fx.header.push_back(line); continue; }
        if (line.size() < 2 || line[1] != ' ') die("malformed line: " + line.substr(0, 60));
        const char c = line[0];
        if (c == 'Q') { if (r || !fx.qLine.empty()) die("Q line must come once, before the records"); fx.qLine = line; continue; }
        if (c == 'K') {
            fx.recs.emplace_back(); r = &fx.recs.back(); r->kLine = line; kv(line, r->k);
            const std::string kind = r->k["kind"];
            if (kind != "C" && kind != "F" && kind != "R") die("unknown record kind '" + kind + "' (" + line.substr(0, 40) + ")");
            for (const char* key : {"kind", "seq", "dim", "pass", "status", "subInfo"})
                if (!r->k.count(key)) die(std::string("K line lacks ") + key + ": " + line.substr(0, 60));
            char* end = 0; r->code = (int)std::strtol(r->k["pass"].c_str(), &end, 10);
            if (r->k["pass"].empty() || *end) die("K line pass= is not an integer: " + line.substr(0, 60));
            r->kind = kind[0];
            continue;
        }
        if (!r) die("record line before the first K line: " + line.substr(0, 40));
        switch (c) {
        case 'G': r->gLine = line; break;
        case 'P': r->pLine = line; kv(line, r->pkv); break;
        case 'E': copy(nums(line, 6), r->eps); r->stateLines.push_back(line); break;
        case 'S': copy(nums(line, 6), r->sig); r->stateLines.push_back(line); break;
        case 'X': copy(nums(line, 6), r->sigEff); r->stateLines.push_back(line); break;
        case 'H': copy(nums(line, 11), r->h); r->stateLines.push_back(line); break;
        case 'I': copy(nums(line, 11), r->im); r->stateLines.push_back(line); break;
        case 'C': copy(nums(line, 2), r->c); r->stateLines.push_back(line); break;
        case 'T': copy(nums(line, 6), r->strain); r->stateLines.push_back(line); break;
        case 'O': copy(nums(line, 17), r->o); asInt(r->o[0], "O status", line); asInt(r->o[1], "O sub", line); r->hasO = true; break;
        case 'D': case 'M': if (r->hasD) die("two tangent lines in one record: " + r->kLine);
                  copy(nums(line, 36), r->D); r->hasD = true; r->dKind = c; break;
        case 'W': copy(nums(line, 12), r->w); r->hasW = true; break;
        case 'Z': { auto v = nums(line, 6); Pass p{asInt(v[0], "Z pass", line), asInt(v[4], "Z status", line), asInt(v[5], "Z sub", line), v[1], v[2], v[3]}; r->z.push_back(p); break; }
        default: die(std::string("unknown line kind '") + c + "' (the fixture is the record-dump subset: Q K G P E S X H I C T O D M W Z)");
        }
    }
    if (fx.recs.empty()) die("the fixture has no records");
    for (const Rec& q : fx.recs) {
        if (q.pLine.empty() || q.stateLines.size() != 7 || !q.hasO || !q.hasW) die("incomplete record: " + q.kLine);
        if (!q.hasD || q.dKind != (q.kind == 'C' ? 'D' : 'M'))
            die(std::string("a ") + q.kind + " record must carry its 6x6 tangent as a " + (q.kind == 'C' ? "D" : "M") + " line: " + q.kLine);
        if (q.kind == 'R' && (q.z.empty() || (q.k.at("dim") != "4" && q.k.at("dim") != "3"))) die("R record without a condensed sigma_22 trace: " + q.kLine);
        if (q.kind != 'R' && !q.z.empty()) die("Z lines on a non-R record: " + q.kLine);
        if (q.k.at("status") != std::to_string((long)q.o[0]) || q.k.at("subInfo") != std::to_string((long)q.o[1]))
            die("K line status=/subInfo= disagree with the O line's status/sub: " + q.kLine);
    }
    return fx;
}

// EVERY Params field from the P line; a key the kernel does not have, or a field the line lacks, is an error.
static void fillParams(const Rec& r, K::Params& p, double& dt) {
    std::map<std::string, std::string> m = r.pkv;
    auto take = [&](const char* key) -> std::string { auto it = m.find(key); if (it == m.end()) die(std::string("P line lacks ") + key + ": " + r.kLine); std::string v = it->second; m.erase(it); return v; };
    auto d = [&](const char* key) { return std::strtod(take(key).c_str(), 0); };
    auto i = [&](const char* key) { return std::atoi(take(key).c_str()); };
    p.E = d("E"); p.nu = d("nu"); p.rho = d("rho"); p.fc = d("fc"); p.ft = d("ft"); p.e = d("e"); p.m0 = d("m0");
    p.Gf = d("Gf"); p.Gc = d("Gc"); p.Df = d("Df"); p.As = d("As"); p.qh0 = d("qh0"); p.Hp = d("Hp");
    p.Ah = d("Ah"); p.Bh = d("Bh"); p.Ch = d("Ch"); p.Dh = d("Dh"); p.epsFc = d("epsFc"); p.eta = d("eta");
    p.implexRmax = d("implexRmax"); p.lch = d("lch"); p.lch_ref = d("lch_ref");
    p.ctTemper = i("ctTemper"); p.tcTemper = i("tcTemper"); p.tensionLaw = i("tensionLaw");
    p.flowPotential = i("flowPotential"); p.maxSubIncr = i("maxSubIncr"); p.compDrive = i("compDrive"); p.implex = i("implex") != 0;
    // the honest-f tolerance is a kernel constant, not a parameter: records written under another one do not replay
    if (d("ftol") != K::F_TOL_HONEST) die("record written with F_TOL_HONEST " + r.pkv.at("ftol") + ", this kernel has another: " + r.kLine);
    dt = d("dt");
    if (!m.empty()) die("P line has a key this testbed does not map: " + m.begin()->first + " (" + r.kLine + ")");
}
static void fillState(const Rec& r, K::State& s) {
    for (int i = 0; i < 6; ++i) { s.eps[i] = r.eps[i]; s.sig[i] = r.sig[i]; s.sigEff[i] = r.sigEff[i]; s.depl[i] = r.im[5 + i]; }
    s.kp = r.h[0]; s.et_max = r.h[1]; s.kdt1 = r.h[2]; s.kdt2 = r.h[3]; s.kdc = r.h[4]; s.kdc1 = r.h[5]; s.kdc2 = r.h[6];
    s.sigtMax = r.h[7]; s.sigcMax = r.h[8]; s.eqc = r.h[9]; s.etPrev = r.h[10];
    s.wt = r.im[0]; s.wc = r.im[1]; s.dwt = r.im[2]; s.dwc = r.im[3]; s.dt_n = r.im[4];
}

struct Ev { int status = 0, sub = 0; double sig[6], se[6], kp = 0, wt = 0, wc = 0, D[6][6]; };
static Ev integrate(const K::Params& mp, const K::State& in, const double s6[6], double dt) {
    Ev e; K::State out;   // LadrunoConcrete3D::integrate: doTangent, hardening on, omega_t/omega_c out
    e.status = K::returnMap(mp, s6, in, out, e.sig, e.se, e.D, true, dt, true, &e.wt, &e.wc);
    e.sub = out.subInfo; e.kp = out.kp; return e;
}
// LadrunoConcrete3D::setTrialStrain, the condensed (PlaneStress / PlateFiber) branch, rule B. Code: 0 sigma_22 not
// converged, 1 last pass failed, 2 both, 3 vanished pivot with the last pass converged; -1 accepted.
struct Loop { int code = -1; std::vector<Pass> z; Ev last; };
static Loop sigma22Loop(const K::Params& mp, const K::State& in, const double strain[6], double eps22, double dt) {
    double s6[6]; for (int i = 0; i < 6; ++i) s6[i] = strain[i]; s6[2] = eps22;
    Loop L; bool lastFailed = false, conv = false, pivot = false; const int maxIt = 30;
    for (int it = 0; it < maxIt; ++it) {
        Ev m = integrate(mp, in, s6, dt); L.last = m; lastFailed = (m.status != 0);
        const double d22 = m.D[2][2];
        double smag = 0.0; for (int i = 0; i < 6; ++i) smag += m.sig[i] * m.sig[i]; smag = std::sqrt(smag);
        L.z.push_back(Pass{it, m.status, m.sub, s6[2], m.sig[2], d22});
        const double tol22 = 1.0e-9 * (smag > 1.0 ? smag : 1.0);
        if (std::fabs(m.sig[2]) <= tol22) { conv = true; break; }
        if (std::fabs(d22) < 1.0e-300) { pivot = true; break; }
        s6[2] -= m.sig[2] / d22;
    }
    if (lastFailed || !conv) L.code = lastFailed ? (conv ? 1 : 2) : (pivot ? 3 : 0);
    return L;
}

static long long ulp(double a, double b) {   // 0 = bit-exact
    if (std::memcmp(&a, &b, 8) == 0) return 0;
    if (!std::isfinite(a) || !std::isfinite(b) || std::signbit(a) != std::signbit(b)) return -1;
    long long ia, ib; std::memcpy(&ia, &a, 8); std::memcpy(&ib, &b, 8); return ia > ib ? ia - ib : ib - ia;
}
static void outcome(const Ev& e, double o[17]) {
    o[0] = e.status; o[1] = e.sub; o[2] = e.kp; o[3] = e.wt; o[4] = e.wc;
    for (int i = 0; i < 6; ++i) { o[5 + i] = e.sig[i]; o[11 + i] = e.se[i]; }
}
static void eig(const double se[6], double w[12]) {
    double A[3][3], ev[3], V[3][3]; K::voigtToMat(se, A); K::eig3sym(A, ev, V);
    for (int i = 0; i < 3; ++i) w[i] = ev[i]; for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) w[3 + 3 * i + j] = V[i][j];
}
static std::string qProbe() {   // the dump's Q line: 21 transcendental probes of the writing binary's libm
    static volatile double qx[3] = {0.37, 1.7, 37.1}, qy[3] = {0.7, -0.9, 2.5}, qz[3] = {0.3, -0.77, 0.999};
    char buf[2048]; int n = std::snprintf(buf, sizeof buf, "Q");
    for (int i = 0; i < 3; i++)
        n += std::snprintf(buf + n, sizeof buf - n, " log=%.17g exp=%.17g cos=%.17g sin=%.17g pow=%.17g atan2=%.17g acos=%.17g",
                           std::log(qx[i]), std::exp(-qx[i]), std::cos(qx[i]), std::sin(qx[i]), std::pow(qx[i], 1.5),
                           std::atan2(qy[i], qz[i]), std::acos(qz[i]));
    return buf;
}

struct Replay { Ev ev; Loop loop; double o[17], w[12]; };
static Replay replay(const Rec& r) {
    K::Params mp; double dt = 0; fillParams(r, mp, dt); K::State in; fillState(r, in);
    Replay x;
    if (r.kind == 'R') { x.loop = sigma22Loop(mp, in, r.strain, r.c[1], dt); x.ev = x.loop.last; }
    else x.ev = integrate(mp, in, r.strain, dt);
    outcome(x.ev, x.o); eig(x.ev.se, x.w);
    return x;
}

static int modeCheck(const Fixture& fx, bool cross) {
    const std::string libm = fx.meta.count("libm") ? fx.meta.at("libm") : "";
    if (libm.empty()) die("the fixture's #@ line names no libm");
    std::printf("fixture: %zu records, recorded with libm=%s kernel-blob=%s; this build: libm=%s\n", fx.recs.size(), libm.c_str(),
                fx.meta.count("kernel-blob") ? fx.meta.at("kernel-blob").c_str() : "?", C3D_REPLAY_LIBM);
    if (libm != C3D_REPLAY_LIBM && !cross)
        die("records written under libm=" + libm + " do not replay bit-exactly under " + C3D_REPLAY_LIBM +
            " (run `record` for record-and-compare on this libm, or --cross to measure the difference)");
    if (cross) std::printf("CROSS-LIBM comparison (diagnostic): a mismatch here is the libm, not the kernel\n");
    int bad = 0; const bool qOk = (qProbe() == fx.qLine);
    std::printf("libm probes (Q line, 21 values): %s\n", qOk ? "bit-identical" : "DIFFER");
    if (!qOk) ++bad;
    std::map<std::string, std::pair<int, int>> byGroup;   // selection class -> (records, bit-exact)
    for (const Rec& r : fx.recs) {
        const Replay x = replay(r);
        std::string why; int nDiff = 0;   // every differing value is counted; the first few are named
        auto note = [&](const std::string& s) { if (nDiff++ < 6) why += " " + s; };
        auto cmp = [&](const std::string& what, double a, double b) {
            const long long u = ulp(a, b);
            if (u) note(what + (u < 0 ? "(sign/nonfinite)" : "(" + std::to_string(u) + " ulp)"));
        };
        auto cmpInt = [&](const std::string& what, long a, long b) { if (a != b) note(what + " " + std::to_string(a) + "!=" + std::to_string(b)); };
        static const char* oName[17] = {"status", "pieces", "kappa", "omega_t", "omega_c", "sig0", "sig1", "sig2", "sig3", "sig4", "sig5", "sigEff0", "sigEff1", "sigEff2", "sigEff3", "sigEff4", "sigEff5"};
        for (int i = 0; i < 17; ++i) { if (i < 2) cmpInt(oName[i], (long)x.o[i], (long)r.o[i]); else cmp(oName[i], x.o[i], r.o[i]); }
        for (int i = 0; i < 36; ++i) cmp(std::string(1, r.dKind) + std::to_string(i / 6) + std::to_string(i % 6), (&x.ev.D[0][0])[i], r.D[i]);
        for (int i = 0; i < 12; ++i) cmp("eig3sym" + std::to_string(i), x.w[i], r.w[i]);
        if (r.kind == 'R') {
            cmpInt("reason", x.loop.code, r.code);
            cmpInt("passes", (long)x.loop.z.size(), (long)r.z.size());
            for (size_t p = 0; p < x.loop.z.size() && p < r.z.size(); ++p) {
                const Pass &a = x.loop.z[p], &b = r.z[p]; const std::string at = "pass" + std::to_string(p) + ".";
                cmpInt(at + "status", a.status, b.status); cmpInt(at + "pieces", a.sub, b.sub);
                cmp(at + "eps22", a.eps22, b.eps22); cmp(at + "s22", a.s22, b.s22); cmp(at + "D22", a.d22, b.d22);
            }
        }
        std::map<std::string, std::string> g; if (!r.gLine.empty()) kv(r.gLine, g);
        auto& cnt = byGroup[g.count("class") ? g["class"] : std::string("?")]; ++cnt.first;
        if (!nDiff) ++cnt.second;
        else { ++bad; std::printf("  MISMATCH %s seq=%s (%s %s): %d values differ:%s%s\n", r.k.at("kind").c_str(), r.k.at("seq").c_str(),
                                  g.count("dump") ? g["dump"].c_str() : "?", g.count("class") ? g["class"].c_str() : "?", nDiff, why.c_str(), nDiff > 6 ? " ..." : ""); }
    }
    for (const auto& kvp : byGroup) std::printf("  %-14s %2d records, bit-exact %2d\n", kvp.first.c_str(), kvp.second.first, kvp.second.second);
    std::printf("REPLAY %s: %d of %zu records bit-exact%s\n", bad ? "FAIL" : "PASS", (int)fx.recs.size() - (bad - (qOk ? 0 : 1)), fx.recs.size(), qOk ? "" : ", libm probes differ");
    return bad ? 1 : 0;
}

static std::string fmt(const double* v, int n) { std::string s; char b[40]; for (int i = 0; i < n; ++i) { std::snprintf(b, sizeof b, " %.17g", v[i]); s += b; } return s; }
static int modeRecord(const Fixture& fx, const char* outPath, const std::string& blob) {
    if (blob.size() != 40 || blob.find_first_not_of("0123456789abcdef") != std::string::npos)
        die("--kernel-blob must be the 40-hex git blob of the kernel header this build compiled (git hash-object)");
    std::ofstream out(outPath); if (!out) die(std::string("cannot write ") + outPath);
    out << "# RE-RECORDED by `concrete3d_replay record` under libm=" << C3D_REPLAY_LIBM << " on kernel-blob=" << blob
        << ": every O/D/M/W/Z line and the K line's outcome fields are THIS build's output; the inputs (P E S X H I C T),\n"
        << "# the selection (G) and the order are the source fixture's. The source fixture's header follows.\n";
    for (const std::string& h : fx.header) {
        if (h.rfind("#@ ", 0) != 0) { out << h << "\n"; continue; }
        out << "#@ libm=" << C3D_REPLAY_LIBM << " kernel-blob=" << blob << " binary=concrete3d_replay-record\n";
    }
    out << qProbe() << "\n";
    for (const Rec& r : fx.recs) {
        const Replay x = replay(r);
        std::map<std::string, std::string> k = r.k;
        if (r.kind == 'R') k["pass"] = std::to_string(x.loop.code);
        k["status"] = std::to_string(x.ev.status); k["subInfo"] = std::to_string(x.ev.sub);   // every kind: they mirror O[0], O[1]
        static const char* order[] = {"kind", "seq", "tag", "ele", "dim", "pass", "status", "subInfo", "cls", "calls", "integ"};
        out << "K"; for (const char* key : order) { auto it = k.find(key); if (it == k.end()) die(std::string("K line lacks ") + key); out << " " << key << "=" << it->second; k.erase(it); }
        if (!k.empty()) die("K line has an unexpected key " + k.begin()->first);
        out << "\n";
        if (!r.gLine.empty()) out << r.gLine << "\n";
        out << r.pLine << "\n"; for (const std::string& s : r.stateLines) out << s << "\n";
        out << "O" << fmt(x.o, 17) << "\n";
        out << (r.kind == 'C' ? "D" : "M") << fmt(&x.ev.D[0][0], 36) << "\n";
        out << "W" << fmt(x.w, 12) << "\n";
        for (const Pass& p : x.loop.z) { const double v[6] = {(double)p.p, p.eps22, p.s22, p.d22, (double)p.status, (double)p.sub}; out << "Z" << fmt(v, 6) << "\n"; }
    }
    std::printf("recorded %zu records under libm=%s to %s\n", fx.recs.size(), C3D_REPLAY_LIBM, outPath);
    return 0;
}

// A glibc tag names the glibc version (linux-glibc-2.39-0ubuntu8.4): the C library this process actually loaded must be
// that version, so a build tagged on one image cannot run as if it were on another.
static void checkRuntimeLibm() {
#ifdef __GLIBC__
    const std::string tag = C3D_REPLAY_LIBM, pre = "linux-glibc-", run = gnu_get_libc_version();
    if (tag.rfind(pre, 0) != 0) die("a glibc build must carry a linux-glibc-<version> tag, not " + tag);
    const std::string ver = tag.substr(pre.size());
    if (ver.compare(0, run.size(), run) != 0 || (ver.size() > run.size() && ver[run.size()] != '-'))
        die("built for glibc " + ver + " but running on glibc " + run);
#endif
}

int main(int argc, char** argv) {
    checkRuntimeLibm();
    if (argc >= 3 && std::strcmp(argv[1], "check") == 0) {
        bool cross = false; for (int i = 3; i < argc; ++i) { if (std::strcmp(argv[i], "--cross") == 0) cross = true; else die(std::string("unknown option ") + argv[i]); }
        return modeCheck(load(argv[2]), cross);
    }
    if (argc == 6 && std::strcmp(argv[1], "record") == 0 && std::strcmp(argv[4], "--kernel-blob") == 0)
        return modeRecord(load(argv[2]), argv[3], argv[5]);
    std::fprintf(stderr, "usage: %s check <fixture> [--cross] | record <fixture> <out> --kernel-blob <sha>\n", argv[0]);
    return 2;
}
