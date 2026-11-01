// concrete3d_replay.cpp -- record-replay testbed for LadrunoConcrete3DKernel.h (apeGmsh #2621, increment F1).
//
// The fixture (concrete3d_replay_records.txt) is a subset of the LADRUNO_C3D_DIAG record dump written by the shipped
// opensees.pyd on the #2602 corbel (fork probe/2621-recdump at ec969aef2, a counters-only instrumentation of
// d14fb6958). A record carries the full committed State, every Params field and the trial strain, and what the binary
// returned. This testbed puts each record back through the kernel OF THIS TREE and requires the recorded answer BIT FOR
// BIT; nothing is compared with a tolerance.
//   C, F  one integrate(): K::returnMap exactly as LadrunoConcrete3D::integrate calls it.
//   R     one refused plane-stress setTrialStrain: the sigma_22 Newton of LadrunoConcrete3D::setTrialStrain (rule B,
//         copied below) from the recorded eps_22 start; the outcome is the last pass's.
// Bit-exactness holds only on the libm the records were written with: Intel's libmmt for the Windows pyd (its Fortran
// objects pull it in by /DEFAULTLIB), glibc for the Linux build. The build defines C3D_REPLAY_LIBM to the libm it
// links (for glibc linux-glibc-<the version ldd --version reports>, checked against the loaded libc at start-up), and
// `check` refuses a fixture recorded under another one, another glibc version included (--cross compares anyway).
// `check` also takes the git blob of the kernel header this build compiled (--kernel-blob, computed by the build
// scripts) and FAILS when the fixture was recorded on another blob. So ANY change to the kernel header fails `check`
// BY DESIGN, and a change that moves results also fails the value compares below (the full returned State included:
// ci.sh proves it on every run with a mutated kernel). Such a PR re-records the fixture on each libm (`record`) and
// commits it; the fixture's diff is the exact list of values the change moved (only the #@ line when none moved), and
// the reviewers judge that diff.
//
//   concrete3d_replay check  <fixture> --kernel-blob <sha> [--cross]  exit 0 iff every record replays bit-exactly
//   concrete3d_replay record <fixture> <out> --kernel-blob <sha>       re-record every OUTPUT on this build (inputs,
//                                                                      order and selection kept), stamped with this
//                                                                      libm and the kernel header's git blob
//   concrete3d_replay inputs-match <a> <b>                             exit 0 iff the two fixtures hold the same
//                                                                      records in the same order with the same inputs
// Exit codes: 0 pass; 1 a compared value (or the kernel blob, or the libm probes) differs; 2 refused at load.
// Builds: concrete3d_replay_ci.sh (g++, glibc) and concrete3d_replay_libmmt.bat (MSVC + Intel libmmt, as the pyd).
//
// FIELD CLASSIFICATION -- every line kind and every field of the fixture. Classes:
//   IN     a kernel input: the replay feeds it to the kernel, and it reaches a compared output
//   OUT    a compared output (the compare named); `record` rewrites it
//   PIN    required at load to equal a dump rule or another field, bit for bit
//   UNREAD an input the kernel does not read in this fixture's modes (implex=0, eta=0, epsFc>0). A coverage limit,
//          stated and ENFORCED: ci.sh's `unread inputs` table perturbs each one by 25 % (0 -> 0.25) and requires the
//          replay to still pass, so a kernel change that starts reading one turns that test red and this table must
//          be updated. (implexRmax and lch_ref are not in it: they are PIN, so a perturbed one is refused at load.)
//   INFO   informational, not compared, with the reason. Across the two fixtures `inputs-match` ties every field
//          marked (x) together, so neither fixture can drop, reorder or edit a record on its own.
// Every number is canonical: a double token must be exactly its own %.17g (what every writer prints), an integer token
// exactly its decimal (no '+', leading zero or -0); anything else is refused at load.
//   #  (header) comment lines: INFO. Only before the first K line; a # line inside the records is refused (it could
//      hide a record or a line). `record` writes its banner as #= lines and drops a source's #= lines.
//   #@ libm=        the libm the outputs were written under: `check` refuses another one (exit 2; --cross waives).
//      kernel-blob= the kernel header's git blob the outputs were written on: compared with --kernel-blob (exit 1).
//      binary=      INFO: who wrote the O/D/W/Z values and the K outcome fields (opensees.pyd@ec969aef2 for the
//                   libmmt fixture; `record` keeps it only when it reproduced every one of those values on the same
//                   libm, otherwise it writes concrete3d_replay-record).
//      M-lines= N-lines=  INFO: who wrote the M and N lines (concrete3d_replay-record).
//   Q  21 libm probes (log exp cos sin pow atan2 acos at three points): OUT, compared with this build's own probes.
//   K  kind=        IN (x): selects the replay (C/F one integrate, R the sigma_22 loop).
//      seq= tag= ele= calls= integ=  INFO (x): the dump's identifiers and provenance counters; the kernel has none.
//      dim=         R: PIN to a condensed element (3 or 4). C/F: INFO (x): the wrapper's element dimension (the
//                   kernel is always 3-D; the record's T is the full 6-component strain).
//      pass=        R: OUT (the rule-B refusal code of the replayed loop). C/F: INFO (x): the sigma_22 pass index of
//                   the recorded integrate(); its strain is T, and pass 0 ties C[1] to T[2] (below).
//      status= subInfo=  PIN to O[0] and O[1] (which are OUT). subInfo also PIN to the dump's domain: -1 when status
//                   != 0, else 0 or (with maxSubIncr > 0) the piece count >= 2; and to N's subInfo (N[40]).
//      cls=         OUT: the dump's return class (0 direct, 1 sub-incremented, 2 floor A, 3 floor B, 4 floor C,
//                   5 direct failure), DERIVED from the replayed returnMap: status 0 -> 0 if subInfo is 0, else 1;
//                   status != 0 -> 4 if maxSubIncr > 0, else 5. This tree's kernel has no class API (the probe's
//                   DiagMap is not in it): the return status and State::subInfo are its whole report, and they
//                   determine every class the probe could write. Class 3 is never assigned by the probe; class 2 needs
//                   a piece f < 2^-maxSubIncr, impossible because every piece and every partial sum is an exact
//                   multiple of 2^-maxSubIncr; a failure with maxSubIncr > 0 leaves only through floor C. Also PIN at
//                   load: the K line's cls must be that formula of its own status, subInfo and the P line's maxSubIncr.
//   G  dump= class= pick=  INFO (x): why the record was picked; domain-checked; class= groups the report.
//   P  E nu fc ft e m0 Gf Df As qh0 Hp Ah Bh Ch Dh epsFc lch ctTemper tcTemper tensionLaw flowPotential compDrive
//      implex eta   IN (x). eta and implex are read as mode switches; eta=0 and implex=0 on every record, so the
//                   Duvaut-Lions and IMPL-EX paths are a coverage limit of this fixture, not exercised.
//      dt           IN (x): read into State::dt_n (N) when dt > 0 (Duvaut-Lions and IMPL-EX are off).
//      maxSubIncr   PIN: 10 iff flowPotential=1, else 0 (the wrapper's rule).
//      implexRmax   PIN: 2 (the dump prints the constant 2.0; the wrapper never sets it). Read only under IMPL-EX.
//      lch_ref      PIN: = lch bit for bit (the dump prints lch twice; integrate sets p.lch_ref = lch). Never read.
//      ftol         PIN: = this kernel's F_TOL_HONEST (a constant, not a parameter).
//      rho          UNREAD (x): the kernel has no use of Params::rho.
//      Gc           UNREAD (x): read only when epsFc <= 0 (epsFcOf); epsFc > 0 on every record.
//   E  committed strain State::eps[6]: IN (x).
//   S  committed nominal stress State::sig[6]: UNREAD (x): returnMap never reads in.sig (it returns from X).
//   X  committed effective stress State::sigEff[6]: IN (x).
//   H  kp et_max kdt1 kdt2 kdc kdc1 kdc2 sigtMax sigcMax eqc etPrev: IN (x) (each but etPrev also reaches its own N field;
//      the kernel sets N.etPrev from this step's drive, and in.etPrev reaches N.eqc and the outputs through the CDPM2
//      compressive history, compDrive = 1 on every record).
//   I  wt wc: IN (x) (into N dwt, dwc).  dwt dwc: UNREAD (x) (IMPL-EX only).  dt_n: IN (x) (into N dt_n when dt <= 0).
//      depl[6]: UNREAD (x) (IMPL-EX only).
//   C  C[0], the committed eps_22: PIN to E[2] (x). C[1], the sigma_22 loop's start: R: IN (x) (the loop starts
//      there). C/F at pass 0: PIN to T[2] (pass 0 runs at the start). C/F at pass > 0: INFO (x): the passes before the
//      recorded one are not in the record, so the path from C[1] to T[2] is not replayed.
//   T  the trial strain: IN (x), except T[2] on an R record: OUT (the eps_22 the loop left in the trial strain).
//   O  status pieces kappa omega_t omega_c sig[6] sigEff[6]: OUT (17 values; status and pieces as integers).
//   D  the 6x6 tangent the pyd dumped (C records): OUT.  M  the same for F and R (last pass), testbed-recorded: OUT.
//   W  eig3sym of the effective stress (3 eigenvalues + the 3x3 vectors), the libm canary: OUT.
//   Z  per pass of an R record: index (PIN to its position), eps22, s22, D22 (OUT), status, pieces (OUT, integers).
//   N  the full State returnMap returned (R: after the last pass), i.e. what the next step would read: eps[6] sig[6]
//      sigEff[6] kp et_max kdt1 kdt2 kdc kdc1 kdc2 sigtMax sigcMax eqc etPrev wt wc dwt dwc dt_n depl[6] subInfo,
//      41 values, every State field: OUT, required on every record. Testbed-recorded under the fixture's libm, like
//      M. On the libmmt fixture they were recorded by a build that reproduced every pyd-written value of every record
//      bit for bit (`record` keeps #@ binary=opensees.pyd@ec969aef2 only then).
// Tangent lines: a C record must carry D, an R or F record M.
// The sigma_22 loop below is a COPY of the wrapper's: the fixture's R records prove the copy matches the wrapper the
// records were written by; a later change to LadrunoConcrete3D::setTrialStrain must be mirrored here.
#include "SRC/material/nD/LadrunoConcrete3DKernel.h"
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
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

static const int N_STATE = 41;
static const char* nName[N_STATE] = {
    "eps0", "eps1", "eps2", "eps3", "eps4", "eps5", "sig0", "sig1", "sig2", "sig3", "sig4", "sig5",
    "sigEff0", "sigEff1", "sigEff2", "sigEff3", "sigEff4", "sigEff5", "kp", "et_max", "kdt1", "kdt2", "kdc", "kdc1",
    "kdc2", "sigtMax", "sigcMax", "eqc", "etPrev", "wt", "wc", "dwt", "dwc", "dt_n",
    "depl0", "depl1", "depl2", "depl3", "depl4", "depl5", "subInfo"};

struct Pass { int p, status, sub; double eps22, s22, d22; };
struct Rec {
    std::string kLine, gLine; std::map<std::string, std::string> k, g;   // K line (verbatim + parsed), G line (selection)
    char kind = 0; int code = 0, status = 0, subInfo = 0, cls = 0;
    std::map<std::string, std::string> pkv; std::string pLine;
    std::vector<std::string> stateLines;                              // E S X H I C T, verbatim
    std::vector<std::string> tTok;                                    // T tokens (inputs-match exempts R's T[2])
    double eps[6], sig[6], sigEff[6], h[11], im[11], c[2], strain[6];
    bool hasO = false, hasD = false, hasW = false, hasN = false; char dKind = 0;   // 'D' (pyd-written) or 'M' (testbed)
    double o[17], D[36], w[12], n[N_STATE];
    std::vector<Pass> z;
    std::set<char> seen;                                              // every line kind but Z comes at most once
};

// A double token must parse whole (no trailing junk, no empty value), must not overflow to +-inf, and must be
// CANONICAL: exactly the %.17g text of the value it parses to (every writer of the fixture prints %.17g), so no two
// spellings of one value (0 / 0.0 / +0 / 1e-3 / 0.001) can pass.
static double num(const std::string& t, const std::string& line) {
    char* end = 0; errno = 0; const double x = std::strtod(t.c_str(), &end);
    if (t.empty() || end == t.c_str() || *end) die("not a number: '" + t + "' in: " + line.substr(0, 60));
    if (errno == ERANGE && std::isinf(x)) die("number overflows a double: " + t + " in: " + line.substr(0, 60));
    char b[40]; std::snprintf(b, sizeof b, "%.17g", x);
    if (t != b) die("not canonical (%.17g gives " + std::string(b) + "): '" + t + "' in: " + line.substr(0, 60));
    return x;
}
// An integer token (status, piece count, pass index, an integer parameter, a K field): exactly the decimal text of an
// int (no '+', no leading zero, no -0, no fraction or exponent), never truncated, within int's range.
static int intTok(const std::string& t, const std::string& what, const std::string& line) {
    char* end = 0; errno = 0; const long long v = std::strtoll(t.c_str(), &end, 10);
    if (t.empty() || end == t.c_str() || *end) die(what + " is not an integer: '" + t + "' in: " + line.substr(0, 60));
    if (errno == ERANGE || v < INT_MIN || v > INT_MAX) die(what + " is out of int range: " + line.substr(0, 60));
    if (std::to_string(v) != t) die(what + " is not a canonical integer: '" + t + "' in: " + line.substr(0, 60));
    return (int)v;
}
// A counter (seq, calls, integ): a canonical non-negative decimal; it may outgrow int.
static void counterTok(const std::string& t, const std::string& what, const std::string& line) {
    if (t.empty() || t.size() > 19 || t.find_first_not_of("0123456789") != std::string::npos || (t.size() > 1 && t[0] == '0'))
        die(what + " is not a canonical whole number: '" + t + "' in: " + line.substr(0, 60));
}
static std::vector<std::string> toks(const std::string& line) {
    std::istringstream is(line.substr(1)); std::vector<std::string> v; std::string t;
    while (is >> t) v.push_back(t);
    return v;
}
// `ints` = the columns that are integers (checked as integer tokens on top of the double check)
static std::vector<double> nums(const std::string& line, size_t want, std::initializer_list<size_t> ints = {}) {
    const std::vector<std::string> t = toks(line); std::vector<double> v;
    for (const std::string& s : t) v.push_back(num(s, line));
    if (v.size() != want) die("line has " + std::to_string(v.size()) + " values, expected " + std::to_string(want) + ": " + line.substr(0, 60));
    for (size_t i : ints) intTok(t[i], std::string(1, line[0]) + "[" + std::to_string(i) + "]", line);
    return v;
}
static void kv(const std::string& line, std::map<std::string, std::string>& m) {
    std::istringstream is(line.substr(1)); std::string t;
    while (is >> t) {
        const size_t p = t.find('=');
        if (p == std::string::npos || p == 0 || p + 1 == t.size()) die("not key=value: " + t + " in: " + line.substr(0, 60));
        if (!m.emplace(t.substr(0, p), t.substr(p + 1)).second) die("duplicate key " + t.substr(0, p) + " in: " + line.substr(0, 60));
    }
}

// The dump's return class, from what this tree's kernel reports (cls= in the classification above).
static int clsOf(int status, int sub, int maxSubIncr) { return status == 0 ? (sub == 0 ? 0 : 1) : (maxSubIncr > 0 ? 4 : 5); }

static void fillParams(const Rec& r, K::Params& p, double& dt);

struct Fixture { std::map<std::string, std::string> meta; std::string qLine; std::vector<Rec> recs; std::vector<std::string> header; };

// needN: every record must carry its N line (check, inputs-match); `record` also reads a source written before N lines.
static Fixture load(const char* path, bool needN = true) {
    std::ifstream in(path); if (!in) die(std::string("cannot open ") + path);
    Fixture fx; std::string line; Rec* r = 0; bool haveMeta = false;
    auto copy = [](const std::vector<double>& v, double* d) { for (size_t i = 0; i < v.size(); ++i) d[i] = v[i]; };
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (r) die("a # line after the first record (it could hide a record or a line): " + line.substr(0, 60));
            if (line.rfind("#@ ", 0) == 0) { if (haveMeta) die("the #@ line must come once, before the records"); haveMeta = true; kv(line.substr(1), fx.meta); }
            fx.header.push_back(line);
            continue;
        }
        if (line.size() < 2 || line[1] != ' ') die("malformed line: " + line.substr(0, 60));
        const char c = line[0];
        if (c == 'Q') { if (r || !fx.qLine.empty()) die("Q line must come once, before the records"); fx.qLine = line; continue; }
        if (c == 'K') {
            fx.recs.emplace_back(); r = &fx.recs.back(); r->kLine = line; kv(line, r->k);
            static const std::set<std::string> kKeys = {"kind", "seq", "tag", "ele", "dim", "pass", "status", "subInfo", "cls", "calls", "integ"};
            for (const std::string& key : kKeys) if (!r->k.count(key)) die("K line lacks " + key + ": " + line.substr(0, 60));
            for (const auto& e : r->k) if (!kKeys.count(e.first)) die("K line has an unexpected key " + e.first + ": " + line.substr(0, 60));
            const std::string kind = r->k.at("kind");
            if (kind != "C" && kind != "F" && kind != "R") die("unknown record kind '" + kind + "' (" + line.substr(0, 40) + ")");
            for (const char* key : {"tag", "ele", "dim"}) intTok(r->k.at(key), std::string("K ") + key, line);
            for (const char* key : {"seq", "calls", "integ"}) counterTok(r->k.at(key), std::string("K ") + key, line);
            r->code = intTok(r->k.at("pass"), "K pass", line);
            r->status = intTok(r->k.at("status"), "K status", line);
            r->subInfo = intTok(r->k.at("subInfo"), "K subInfo", line);
            r->cls = intTok(r->k.at("cls"), "K cls", line);
            r->kind = kind[0];
            continue;
        }
        if (!r) die("record line before the first K line: " + line.substr(0, 40));
        if (c != 'Z' && !r->seen.insert(c == 'M' ? 'D' : c).second) die(std::string("two ") + c + " lines in one record: " + r->kLine);
        switch (c) {
        case 'G': {
            r->gLine = line; kv(line, r->g);
            if (r->g.size() != 3 || !r->g.count("dump") || !r->g.count("class") || !r->g.count("pick")) die("G line must carry exactly dump= class= pick=: " + line);
            const std::string& d = r->g.at("dump"); const std::string& cl = r->g.at("class");
            if (d != "unset" && d != "COMPATIBLE") die("G dump= is not unset|COMPATIBLE: " + line);
            if (cl != "refusal" && cl != "floor-C" && cl != "sub-increment" && cl != "kink" && cl != "direct") die("G class= is not a record class: " + line);
            break; }
        case 'P': r->pLine = line; kv(line, r->pkv); break;
        case 'E': copy(nums(line, 6), r->eps); r->stateLines.push_back(line); break;
        case 'S': copy(nums(line, 6), r->sig); r->stateLines.push_back(line); break;
        case 'X': copy(nums(line, 6), r->sigEff); r->stateLines.push_back(line); break;
        case 'H': copy(nums(line, 11), r->h); r->stateLines.push_back(line); break;
        case 'I': copy(nums(line, 11), r->im); r->stateLines.push_back(line); break;
        case 'C': copy(nums(line, 2), r->c); r->stateLines.push_back(line); break;
        case 'T': copy(nums(line, 6), r->strain); r->tTok = toks(line); r->stateLines.push_back(line); break;
        case 'O': copy(nums(line, 17, {0, 1}), r->o); r->hasO = true; break;
        case 'D': case 'M': copy(nums(line, 36), r->D); r->hasD = true; r->dKind = c; break;
        case 'W': copy(nums(line, 12), r->w); r->hasW = true; break;
        case 'N': copy(nums(line, N_STATE, {N_STATE - 1}), r->n); r->hasN = true; break;
        case 'Z': { auto v = nums(line, 6, {0, 4, 5}); Pass p{(int)v[0], (int)v[4], (int)v[5], v[1], v[2], v[3]};
                    if (p.p != (int)r->z.size()) die("Z pass index " + std::to_string(p.p) + " is not its position " + std::to_string(r->z.size()) + ": " + r->kLine);
                    r->z.push_back(p); break; }
        default: die(std::string("unknown line kind '") + c + "' (the fixture is the record-dump subset: Q K G P E S X H I C T O D M W Z N)");
        }
    }
    if (fx.recs.empty()) die("the fixture has no records");
    if (fx.qLine.empty()) die("the fixture has no Q line (libm probes)");
    for (const char* key : {"libm", "kernel-blob"}) if (!fx.meta.count(key)) die(std::string("the #@ line lacks ") + key);
    for (const auto& e : fx.meta)
        if (e.first != "libm" && e.first != "kernel-blob" && e.first != "binary" && e.first != "M-lines" && e.first != "N-lines")
            die("the #@ line has an unexpected key " + e.first);
    { const std::string& b = fx.meta.at("kernel-blob");
      if (b.size() != 40 || b.find_first_not_of("0123456789abcdef") != std::string::npos) die("the #@ kernel-blob is not a 40-hex git blob: " + b); }
    for (const Rec& q : fx.recs) {
        for (char c : std::string(needN ? "GPESXHICTOWN" : "GPESXHICTOW")) if (!q.seen.count(c)) die(std::string("record lacks its ") + c + " line: " + q.kLine);
        if (std::memcmp(&q.c[0], &q.eps[2], sizeof(double)) != 0) die("C line's committed eps_22 is not E[2] bit for bit: " + q.kLine);
        if (q.kind != 'R' && q.code == 0 && std::memcmp(&q.c[1], &q.strain[2], sizeof(double)) != 0)
            die("a pass-0 C/F record's loop start C[1] is not its trial T[2] bit for bit: " + q.kLine);
        K::Params pp; double dt = 0; fillParams(q, pp, dt);   // every P line is validated before anything replays
        if (!q.hasD || q.dKind != (q.kind == 'C' ? 'D' : 'M'))
            die(std::string("a ") + q.kind + " record must carry its 6x6 tangent as a " + (q.kind == 'C' ? "D" : "M") + " line: " + q.kLine);
        if (q.kind == 'R' && (q.z.empty() || (q.k.at("dim") != "4" && q.k.at("dim") != "3"))) die("R record without a condensed sigma_22 trace: " + q.kLine);
        if (q.kind != 'R' && !q.z.empty()) die("Z lines on a non-R record: " + q.kLine);
        if (q.status != (int)q.o[0] || q.subInfo != (int)q.o[1])
            die("K line status=/subInfo= disagree with the O line's status/sub: " + q.kLine);
        const bool subOk = q.status != 0 ? q.subInfo == -1 : (q.subInfo == 0 || (q.subInfo >= 2 && pp.maxSubIncr > 0));
        if (!subOk || q.cls != clsOf(q.status, q.subInfo, pp.maxSubIncr))
            die("K cls=" + std::to_string(q.cls) + " subInfo=" + std::to_string(q.subInfo) + " is not what the dump writes for status=" +
                std::to_string(q.status) + " (maxSubIncr=" + std::to_string(pp.maxSubIncr) + "): " + q.kLine);
        if (q.hasN && (int)q.n[N_STATE - 1] != q.subInfo) die("N line's subInfo is not the K line's subInfo: " + q.kLine);
    }
    return fx;
}

// EVERY Params field from the P line; a key the kernel does not have, or a field the line lacks, is an error.
static void fillParams(const Rec& r, K::Params& p, double& dt) {
    std::map<std::string, std::string> m = r.pkv;
    auto take = [&](const char* key) -> std::string { auto it = m.find(key); if (it == m.end()) die(std::string("P line lacks ") + key + ": " + r.kLine); std::string v = it->second; m.erase(it); return v; };
    auto d = [&](const char* key) { return num(take(key), r.pLine); };
    auto i = [&](const char* key) { return intTok(take(key), std::string("P ") + key, r.pLine); };
    p.E = d("E"); p.nu = d("nu"); p.rho = d("rho"); p.fc = d("fc"); p.ft = d("ft"); p.e = d("e"); p.m0 = d("m0");
    p.Gf = d("Gf"); p.Gc = d("Gc"); p.Df = d("Df"); p.As = d("As"); p.qh0 = d("qh0"); p.Hp = d("Hp");
    p.Ah = d("Ah"); p.Bh = d("Bh"); p.Ch = d("Ch"); p.Dh = d("Dh"); p.epsFc = d("epsFc"); p.eta = d("eta");
    p.implexRmax = d("implexRmax"); p.lch = d("lch"); p.lch_ref = d("lch_ref");
    p.ctTemper = i("ctTemper"); p.tcTemper = i("tcTemper"); p.tensionLaw = i("tensionLaw");
    p.flowPotential = i("flowPotential"); p.maxSubIncr = i("maxSubIncr"); p.compDrive = i("compDrive");
    const int implex = i("implex"); p.implex = implex != 0;
    // only values LadrunoConcrete3D's option parser can produce (-ctTemper none|alphat|proj, -tcTemper none|proj,
    // -tensionLaw exp|bilinear, -flowPotential legacy|cdpm2, -compressionDrive legacy|cdpm2, -implex, maxSubIncr = 10
    // iff the CDPM2 potential): anything else would replay as some other mode without a word
    auto in = [&](const char* key, int v, std::initializer_list<int> ok) {
        for (int o : ok) if (v == o) return;
        die(std::string("P ") + key + "=" + std::to_string(v) + " is not a value the wrapper writes: " + r.kLine);
    };
    in("ctTemper", p.ctTemper, {0, 1, 2}); in("tcTemper", p.tcTemper, {0, 2}); in("tensionLaw", p.tensionLaw, {0, 1});
    in("flowPotential", p.flowPotential, {0, 1}); in("compDrive", p.compDrive, {0, 1}); in("implex", implex, {0, 1});
    in("maxSubIncr", p.maxSubIncr, {p.flowPotential == 1 ? 10 : 0});
    // the dump's two constant rules: LadrunoConcrete3D::diagRecord (ec969aef2) prints the literal 2.0 for implexRmax
    // and lch for both lch and lch_ref (integrate sets p.lch_ref = lch and never sets implexRmax)
    { const double two = 2.0; if (std::memcmp(&p.implexRmax, &two, sizeof(double)) != 0) die("P implexRmax=" + r.pkv.at("implexRmax") + " is not the dump's constant 2: " + r.kLine); }
    if (std::memcmp(&p.lch_ref, &p.lch, sizeof(double)) != 0) die("P lch_ref is not lch bit for bit (the dump prints lch twice): " + r.kLine);
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
// every K::State field, in the N line's order (the input lines' order: E S X H I, then subInfo)
static void stateOut(const K::State& s, double n[N_STATE]) {
    for (int i = 0; i < 6; ++i) { n[i] = s.eps[i]; n[6 + i] = s.sig[i]; n[12 + i] = s.sigEff[i]; n[34 + i] = s.depl[i]; }
    n[18] = s.kp; n[19] = s.et_max; n[20] = s.kdt1; n[21] = s.kdt2; n[22] = s.kdc; n[23] = s.kdc1; n[24] = s.kdc2;
    n[25] = s.sigtMax; n[26] = s.sigcMax; n[27] = s.eqc; n[28] = s.etPrev;
    n[29] = s.wt; n[30] = s.wc; n[31] = s.dwt; n[32] = s.dwc; n[33] = s.dt_n; n[40] = s.subInfo;
}

struct Ev { int status = 0, sub = 0; double sig[6], se[6], kp = 0, wt = 0, wc = 0, D[6][6]; K::State out; };
static Ev integrate(const K::Params& mp, const K::State& in, const double s6[6], double dt) {
    Ev e;   // LadrunoConcrete3D::integrate: doTangent, hardening on, omega_t/omega_c out; e.out = the returned State
    e.status = K::returnMap(mp, s6, in, e.out, e.sig, e.se, e.D, true, dt, true, &e.wt, &e.wc);
    e.sub = e.out.subInfo; e.kp = e.out.kp; return e;
}
// LadrunoConcrete3D::setTrialStrain, the condensed (PlaneStress / PlateFiber) branch, rule B. Code: 0 sigma_22 not
// converged, 1 last pass failed, 2 both, 3 vanished pivot with the last pass converged; -1 accepted.
struct Loop { int code = -1; std::vector<Pass> z; Ev last; double eps22End = 0; };   // eps22End: s6[2] at exit
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
    L.eps22End = s6[2];
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

struct Replay { Ev ev; Loop loop; double o[17], w[12], n[N_STATE]; int cls = 0; };
static Replay replay(const Rec& r) {
    K::Params mp; double dt = 0; fillParams(r, mp, dt); K::State in; fillState(r, in);
    Replay x;
    if (r.kind == 'R') { x.loop = sigma22Loop(mp, in, r.strain, r.c[1], dt); x.ev = x.loop.last; }
    else x.ev = integrate(mp, in, r.strain, dt);
    outcome(x.ev, x.o); eig(x.ev.se, x.w); stateOut(x.ev.out, x.n);
    x.cls = clsOf(x.ev.status, x.ev.sub, mp.maxSubIncr);
    return x;
}

// Every compared value of one record (withN: the N line too); returns how many differ, the first few named in `why`.
static int compare(const Rec& r, const Replay& x, std::string& why, bool withN = true) {
    int nDiff = 0;
    auto note = [&](const std::string& s) { if (nDiff++ < 6) why += " " + s; };
    auto cmp = [&](const std::string& what, double a, double b) {
        const long long u = ulp(a, b);
        if (u) note(what + (u < 0 ? "(sign/nonfinite)" : "(" + std::to_string(u) + " ulp)"));
    };
    auto cmpInt = [&](const std::string& what, long a, long b) { if (a != b) note(what + " " + std::to_string(a) + "!=" + std::to_string(b)); };
    static const char* oName[17] = {"status", "pieces", "kappa", "omega_t", "omega_c", "sig0", "sig1", "sig2", "sig3", "sig4", "sig5", "sigEff0", "sigEff1", "sigEff2", "sigEff3", "sigEff4", "sigEff5"};
    for (int i = 0; i < 17; ++i) { if (i < 2) cmpInt(oName[i], (long)x.o[i], (long)r.o[i]); else cmp(oName[i], x.o[i], r.o[i]); }
    cmpInt("cls", x.cls, r.cls);
    for (int i = 0; i < 36; ++i) cmp(std::string(1, r.dKind) + std::to_string(i / 6) + std::to_string(i % 6), (&x.ev.D[0][0])[i], r.D[i]);
    for (int i = 0; i < 12; ++i) cmp("eig3sym" + std::to_string(i), x.w[i], r.w[i]);
    if (withN) for (int i = 0; i < N_STATE; ++i) {
        if (i == N_STATE - 1) cmpInt(std::string("N.") + nName[i], (long)x.n[i], (long)r.n[i]);
        else cmp(std::string("N.") + nName[i], x.n[i], r.n[i]);
    }
    if (r.kind == 'R') {
        cmpInt("reason", x.loop.code, r.code);
        cmpInt("passes", (long)x.loop.z.size(), (long)r.z.size());
        cmp("T2end", x.loop.eps22End, r.strain[2]);
        for (size_t p = 0; p < x.loop.z.size() && p < r.z.size(); ++p) {
            const Pass &a = x.loop.z[p], &b = r.z[p]; const std::string at = "pass" + std::to_string(p) + ".";
            cmpInt(at + "status", a.status, b.status); cmpInt(at + "pieces", a.sub, b.sub);
            cmp(at + "eps22", a.eps22, b.eps22); cmp(at + "s22", a.s22, b.s22); cmp(at + "D22", a.d22, b.d22);
        }
    }
    return nDiff;
}

static void checkBlobArg(const std::string& blob) {
    if (blob.size() != 40 || blob.find_first_not_of("0123456789abcdef") != std::string::npos)
        die("--kernel-blob must be the 40-hex git blob of the kernel header this build compiled (git hash-object)");
}

static int modeCheck(const Fixture& fx, bool cross, const std::string& blob) {
    const std::string& libm = fx.meta.at("libm");
    std::printf("fixture: %zu records, recorded with libm=%s kernel-blob=%s; this build: libm=%s kernel-blob=%s\n", fx.recs.size(),
                libm.c_str(), fx.meta.at("kernel-blob").c_str(), C3D_REPLAY_LIBM, blob.c_str());
    if (libm != C3D_REPLAY_LIBM && !cross)
        die("records written under libm=" + libm + " do not replay bit-exactly under " + C3D_REPLAY_LIBM +
            " (run `record` for record-and-compare on this libm, or --cross to measure the difference)");
    if (cross) std::printf("CROSS-LIBM comparison (diagnostic): a mismatch here is the libm, not the kernel\n");
    const bool blobOk = (fx.meta.at("kernel-blob") == blob);
    if (!blobOk) std::printf("KERNEL-BLOB MISMATCH: the fixture was recorded on kernel-blob %s, this build compiled %s: re-record (`record`)\n",
                             fx.meta.at("kernel-blob").c_str(), blob.c_str());
    int bad = 0; const bool qOk = (qProbe() == fx.qLine);
    std::printf("libm probes (Q line, 21 values): %s\n", qOk ? "bit-identical" : "DIFFER");
    std::map<std::string, std::pair<int, int>> byGroup;   // selection class -> (records, bit-exact)
    for (const Rec& r : fx.recs) {
        std::string why; const int nDiff = compare(r, replay(r), why);
        const std::map<std::string, std::string>& g = r.g;
        auto& cnt = byGroup[g.at("class")]; ++cnt.first;
        if (!nDiff) ++cnt.second;
        else { ++bad; std::printf("  MISMATCH %s seq=%s (%s %s): %d values differ:%s%s\n", r.k.at("kind").c_str(), r.k.at("seq").c_str(),
                                  g.at("dump").c_str(), g.at("class").c_str(), nDiff, why.c_str(), nDiff > 6 ? " ..." : ""); }
    }
    for (const auto& kvp : byGroup) std::printf("  %-14s %2d records, bit-exact %2d\n", kvp.first.c_str(), kvp.second.first, kvp.second.second);
    const bool fail = bad || !qOk || !blobOk;
    std::printf("REPLAY %s: %d of %zu records bit-exact%s%s\n", fail ? "FAIL" : "PASS", (int)fx.recs.size() - bad, fx.recs.size(),
                qOk ? "" : ", libm probes differ", blobOk ? "" : ", kernel-blob differs");
    return fail ? 1 : 0;
}

static std::string fmt(const double* v, int n) { std::string s; char b[40]; for (int i = 0; i < n; ++i) { std::snprintf(b, sizeof b, " %.17g", v[i]); s += b; } return s; }
static int modeRecord(const Fixture& fx, const char* outPath, const std::string& blob) {
    // #@ binary= stays the source's only when this build reproduces every source output on the source's libm (an N line
    // the source lacks is not a source output)
    bool reproduced = (fx.meta.at("libm") == C3D_REPLAY_LIBM);
    std::vector<Replay> xs;
    for (const Rec& r : fx.recs) {
        xs.push_back(replay(r));
        std::string why;
        if (reproduced && compare(r, xs.back(), why, r.hasN)) reproduced = false;
    }
    const bool keep = reproduced && fx.meta.count("binary");
    const std::string binary = keep ? fx.meta.at("binary") : "concrete3d_replay-record";
    std::ofstream out(outPath, std::ios::binary); if (!out) die(std::string("cannot write ") + outPath);   // LF on every OS
    out << "#= RE-RECORDED by `concrete3d_replay record` under libm=" << C3D_REPLAY_LIBM << " on kernel-blob=" << blob << ".\n"
        << "#= OUTPUTS, this build's: every O, D/M, W, Z and N line, the K line's status=, subInfo= and cls=, an R record's\n"
        << "#= pass= and T[2]. INPUTS and selection, the source fixture's, verbatim: G P E S X H I C, T (but an R record's\n"
        << "#= T[2]), the K line's kind= seq= tag= ele= dim= calls= integ= and a C/F record's pass=, and the record order.\n"
        << "#= Every field's class (input, compared output, pinned, unread, informational) is the FIELD CLASSIFICATION table\n"
        << "#= at the top of tests/_testbed/concrete3d_replay.cpp. #@ binary= is "
        << (keep ? "the source's: this build reproduced every source\n#= output bit for bit on the source's libm."
                 : "this build: it did not reproduce every\n#= source output on the source's libm (or the source names no binary).")
        << " The source fixture's header follows.\n";
    for (const std::string& h : fx.header) {
        if (h.rfind("#=", 0) == 0) continue;   // an earlier `record` banner
        if (h.rfind("#@ ", 0) != 0) { out << h << "\n"; continue; }
        out << "#@ libm=" << C3D_REPLAY_LIBM << " kernel-blob=" << blob << " binary=" << binary
            << " M-lines=concrete3d_replay-record N-lines=concrete3d_replay-record\n";
    }
    out << qProbe() << "\n";
    for (size_t ri = 0; ri < fx.recs.size(); ++ri) {
        const Rec& r = fx.recs[ri]; const Replay& x = xs[ri];
        std::map<std::string, std::string> k = r.k;
        if (r.kind == 'R') k["pass"] = std::to_string(x.loop.code);
        k["status"] = std::to_string(x.ev.status); k["subInfo"] = std::to_string(x.ev.sub);   // every kind: they mirror O[0], O[1]
        k["cls"] = std::to_string(x.cls);
        static const char* order[] = {"kind", "seq", "tag", "ele", "dim", "pass", "status", "subInfo", "cls", "calls", "integ"};
        out << "K"; for (const char* key : order) { auto it = k.find(key); if (it == k.end()) die(std::string("K line lacks ") + key); out << " " << key << "=" << it->second; k.erase(it); }
        if (!k.empty()) die("K line has an unexpected key " + k.begin()->first);
        out << "\n";
        out << r.gLine << "\n";
        out << r.pLine << "\n";
        for (const std::string& s : r.stateLines) {
            if (r.kind == 'R' && s[0] == 'T') { double t[6]; for (int i = 0; i < 6; ++i) t[i] = r.strain[i]; t[2] = x.loop.eps22End; out << "T" << fmt(t, 6) << "\n"; }
            else out << s << "\n";
        }
        out << "O" << fmt(x.o, 17) << "\n";
        out << (r.kind == 'C' ? "D" : "M") << fmt(&x.ev.D[0][0], 36) << "\n";
        out << "W" << fmt(x.w, 12) << "\n";
        for (const Pass& p : x.loop.z) { const double v[6] = {(double)p.p, p.eps22, p.s22, p.d22, (double)p.status, (double)p.sub}; out << "Z" << fmt(v, 6) << "\n"; }
        out << "N" << fmt(x.n, N_STATE) << "\n";
    }
    std::printf("recorded %zu records under libm=%s to %s (binary=%s)\n", fx.recs.size(), C3D_REPLAY_LIBM, outPath, binary.c_str());
    return 0;
}

// The two fixtures (libmmt and glibc) must hold the same records, in the same order, with the same inputs: G P E S X H
// I C verbatim, T token for token (an R record's T[2] is an output), and the K keys that are not outputs.
static int modeInputsMatch(const Fixture& a, const Fixture& b) {
    int bad = 0;
    auto diff = [&](size_t i, const std::string& what) { if (bad++ < 10) std::printf("  INPUTS DIFFER at record %zu (seq=%s): %s\n", i, a.recs[i].k.at("seq").c_str(), what.c_str()); };
    if (a.recs.size() != b.recs.size()) { std::printf("INPUTS DIFFER: %zu records vs %zu\n", a.recs.size(), b.recs.size()); return 1; }
    for (size_t i = 0; i < a.recs.size(); ++i) {
        const Rec &p = a.recs[i], &q = b.recs[i];
        for (const char* key : {"kind", "seq", "tag", "ele", "dim", "calls", "integ"})
            if (p.k.at(key) != q.k.at(key)) diff(i, std::string("K ") + key + "=" + p.k.at(key) + " vs " + q.k.at(key));
        if (p.kind != 'R' && p.k.at("pass") != q.k.at("pass")) diff(i, "K pass=" + p.k.at("pass") + " vs " + q.k.at("pass"));
        if (p.gLine != q.gLine) diff(i, "G line");
        if (p.pLine != q.pLine) diff(i, "P line");
        for (size_t s = 0; s < p.stateLines.size(); ++s) {   // E S X H I C T, each exactly once (load)
            if (p.stateLines[s][0] != q.stateLines[s][0]) { diff(i, "state line order"); break; }
            if (p.stateLines[s][0] != 'T' && p.stateLines[s] != q.stateLines[s]) diff(i, std::string(1, p.stateLines[s][0]) + " line");
        }
        for (size_t t = 0; t < 6; ++t) if (!(p.kind == 'R' && t == 2) && p.tTok[t] != q.tTok[t]) diff(i, "T[" + std::to_string(t) + "]");
    }
    std::printf("INPUTS %s: %zu records%s\n", bad ? "DIFFER" : "MATCH", a.recs.size(), bad ? "" : ", same order, same inputs");
    return bad ? 1 : 0;
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
        bool cross = false; std::string blob;
        for (int i = 3; i < argc; ++i) {
            if (std::strcmp(argv[i], "--cross") == 0) cross = true;
            else if (std::strcmp(argv[i], "--kernel-blob") == 0 && i + 1 < argc) blob = argv[++i];
            else die(std::string("unknown option ") + argv[i]);
        }
        if (blob.empty()) die("check needs --kernel-blob <the git blob of the kernel header this build compiled>");
        checkBlobArg(blob);
        return modeCheck(load(argv[2]), cross, blob);
    }
    if (argc == 6 && std::strcmp(argv[1], "record") == 0 && std::strcmp(argv[4], "--kernel-blob") == 0) {
        checkBlobArg(argv[5]);
        return modeRecord(load(argv[2], false), argv[3], argv[5]);
    }
    if (argc == 4 && std::strcmp(argv[1], "inputs-match") == 0) return modeInputsMatch(load(argv[2]), load(argv[3]));
    std::fprintf(stderr, "usage: %s check <fixture> --kernel-blob <sha> [--cross] | record <fixture> <out> --kernel-blob <sha> | inputs-match <a> <b>\n", argv[0]);
    return 2;
}
