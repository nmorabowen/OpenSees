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
//   concrete3d_replay fields [--header|--rows|--check <this .cpp>]     the field table (below); --rows its generated
//                                                                      rows; --check exit 0 iff the header block is it
//   concrete3d_replay contract <fixture> --kernel-blob <sha>           exit 0 iff every generated row gives its
//                                                                      class's outcome (the two-way contract below)
// Exit codes: 0 pass; 1 a compared value (or the kernel blob, or the libm probes) differs; 2 refused at load.
// Builds: concrete3d_replay_ci.sh (g++, glibc) and concrete3d_replay_libmmt.bat (MSVC + Intel libmmt, as the pyd).
//
// FIELD CLASSIFICATION -- every line kind and every field of the fixture. Classes:
//   IN     a kernel input: the replay feeds it to the kernel, and it reaches a compared output
//   OUT    a compared output (the compare named); `record` rewrites it
//   PIN    required at load to equal a dump rule or another field, bit for bit
//   UNREAD an input the kernel does not read in this fixture's modes (implex=0, eta=0, epsFc>0). A coverage limit,
//          stated and ENFORCED (below)
//   INFO   informational, not compared, with the reason.
// Across the two fixtures `inputs-match` ties the G P E S X H I C T lines (but an R record's T[2]) and the K keys kind
// seq tag ele dim calls integ (and a C/F record's pass) together, so neither fixture can drop, reorder or edit a record
// on its own.
// Every number is canonical: a double token must be exactly its own %.17g (what every writer prints), an integer token
// exactly its decimal (no '+', leading zero or -0); anything else is refused at load.
// The table below is GENERATED: it is the output of `concrete3d_replay fields`, which prints FIELDS (the field table
// in this file, its single source), and `fields --check <this file>` fails when the two differ (ci.sh and the .bat run
// it). The loader takes every numeric line's column count, every key set and every domain from FIELDS, so a field
// cannot be added unclassified. It is a TWO-WAY CONTRACT: `concrete3d_replay contract <fixture>` generates one row per
// column of every entry (per domain value for {}), perturbs that field on every record in its scope (with its ties),
// and requires the class's outcome: an IN row must FAIL (exit 1, its first mismatch naming a compared field) or be
// refused at load where the entry says so; an OUT row must fail naming its own compare; a PIN row must be refused; an
// UNREAD or INFO row must still PASS on every record. So a later edit that stops reading or comparing an IN field turns its row
// red, and a kernel that starts reading an UNREAD one turns that row red. The rows run in-process through the same
// load / replay / compare as `check` (a refusal is exit 2, a difference exit 1); a record whose input lines a row did
// not touch reuses the base fixture's replay (the replay is a pure function of them), and a row that must fail stops
// at its first mismatching record. ci.sh --contract runs them on the glibc fixture (its own workflow step), the .bat's
// `contract` on the libmmt fixture.
//   ---- BEGIN FIELD TABLE (the output of `concrete3d_replay fields`; `fields --check` compares) ----
//   ln cols          class  scope pert    tie            expect -- note
//   #  text          INFO   all   -       -              -- header comments: only before the first K line
//   #@ libm          PIN    all   =other  -              refuse:do not replay bit-exactly under -- the outputs' libm; --cross waives
//   #@ kernel-blob   OUT    all   ~last   -              fail:KERNEL-BLOB MISMATCH -- the kernel header blob, vs --kernel-blob
//   #@ binary        INFO   all   =x      -              -- who wrote the O/D/W/Z values and the K outcome fields (note 2)
//   #@ M-lines       INFO   all   =x      -              -- who wrote the M lines
//   #@ N-lines       INFO   all   =x      -              -- who wrote the N lines
//   Q  0-20          OUT    all   x1.25   -              fail:libm probes (Q line, 21 values): DIFFER -- 21 libm probes, vs this build's
//   K  kind          IN     R     =F      -              refuse:Z lines on a non-R record -- C/F: one integrate, R: the sigma_22 loop; domain C,F,R
//   K  seq           INFO   all   +1      -              -- the dump's identifiers and counters: the kernel has none
//   K  tag           INFO   all   +1      -              -- (as seq)
//   K  ele           INFO   all   +1      -              -- (as seq)
//   K  calls         INFO   all   +1      -              -- (as seq)
//   K  integ         INFO   all   +1      -              -- (as seq)
//   K  dim           INFO   CF    +1      -              -- the wrapper's element dimension: the kernel is always 3-D
//   K  dim           PIN    R     =5      -              refuse:without a condensed sigma_22 trace -- a condensed element: 3 or 4
//   K  pass          OUT    R     +1      -              fail: reason  -- the rule-B refusal code of the replayed loop
//   K  pass          INFO   CF    +1      -              -- the sigma_22 pass of the recorded integrate(); pass 0 ties C[1] to T[2]
//   K  status        PIN    all   +1      -              refuse:disagree with the O line -- = O[0]
//   K  subInfo       PIN    all   +1      -              refuse:disagree with the O line -- = O[1] and N[40], in its domain (note 1)
//   K  cls           OUT    all   +1      -              refuse:is not what the dump writes -- derived (note 1); its load pin refuses first
//   G  dump          INFO   all   =unset  -              -- why the record was picked; domain unset,COMPATIBLE
//   G  class         INFO   all   =direct -              -- groups the report; domain refusal,floor-C,sub-increment,kink,direct
//   G  pick          INFO   all   =x      -              -- the selection rule
//   P  E             IN     all   x1.25   -
//   P  nu            IN     all   x1.25   -
//   P  rho           UNREAD all   x1.25   -              -- the kernel has no use of Params::rho
//   P  fc            IN     all   x1.25   -
//   P  ft            IN     all   x1.25   -
//   P  e             IN     all   x1.25   -
//   P  m0            IN     all   x1.25   -
//   P  Gf            IN     all   x1.25   -
//   P  Gc            UNREAD all   x1.25   -              -- read only when epsFc <= 0 (epsFcOf); epsFc > 0 on every record
//   P  Df            IN     all   x1.25   -
//   P  As            IN     all   x1.25   -
//   P  qh0           IN     all   x1.25   -
//   P  Hp            IN     all   x1.25   -
//   P  Ah            IN     all   x1.25   -
//   P  Bh            IN     all   x1.25   -
//   P  Ch            IN     all   x1.25   -
//   P  Dh            IN     all   x1.25   -
//   P  ctTemper      IN     all   {}      -              -- -ctTemper none|alphat|proj; domain 0,1,2
//   P  tcTemper      IN     all   {}      -              -- -tcTemper none|proj; domain 0,2
//   P  tensionLaw    IN     all   {}      -              -- -tensionLaw exp|bilinear; domain 0,1
//   P  epsFc         IN     all   x1.25   -
//   P  flowPotential IN     cls0  =0      P.maxSubIncr=0 -- cls 0 only: elsewhere the cls pin refuses it; domain 0,1
//   P  maxSubIncr    PIN    all   =0      -              refuse:is not a value the wrapper writes -- 10 iff flowPotential=1, else 0
//   P  compDrive     IN     all   {}      -              -- -compressionDrive legacy|cdpm2; domain 0,1
//   P  eta           IN     all   x1.25   -              -- 0 on every record: 0 -> 0.25 switches Duvaut-Lions on
//   P  implex        IN     all   {}      -              -- 0 on every record: the IMPL-EX path is a coverage limit; domain 0,1
//   P  implexRmax    PIN    all   x1.25   -              refuse:is not the dump's constant 2 -- the dump's constant; IMPL-EX only
//   P  lch           IN     all   x1.25   P.lch_ref
//   P  lch_ref       PIN    all   x1.25   -              refuse:is not lch bit for bit -- the dump prints lch twice; never read
//   P  ftol          PIN    all   x1.25   -              refuse:this kernel has another -- = this kernel's F_TOL_HONEST (a constant)
//   P  dt            IN     dt>0  x1.25   -              -- into N.dt_n
//   P  dt            UNREAD dt<0  x1.25   -              -- only its sign is read (Duvaut-Lions and IMPL-EX off)
//   E  0-1           IN     all   x1.25   -              -- committed strain State::eps
//   E  2             IN     all   x1.25   C.0
//   E  3-5           IN     all   x1.25   -
//   S  0-5           UNREAD all   x1.25   -              -- State::sig: returnMap never reads in.sig (it returns from X)
//   X  0-5           IN     all   x1.25   -              -- committed effective stress State::sigEff
//   H  0-10          IN     all   x1.25   -              -- kp et_max kdt1 kdt2 kdc kdc1 kdc2 sigtMax sigcMax eqc etPrev
//   I  0-1           IN     all   x1.25   -              -- wt wc, into N.dwt N.dwc
//   I  2-3           UNREAD all   x1.25   -              -- dwt dwc: IMPL-EX only
//   I  4             IN     dt<=0 x1.25   -              -- dt_n, into N.dt_n when dt <= 0
//   I  4             UNREAD dt>0  x1.25   -              -- dt_n: N.dt_n = dt when dt > 0
//   I  5-10          UNREAD all   x1.25   -              -- depl[6]: IMPL-EX only
//   C  0             PIN    all   x1.25   -              refuse:is not E[2] bit for bit -- the committed eps_22 = E[2]
//   C  1             IN     R     x1.25   -              -- the sigma_22 loop's start
//   C  1             PIN    CF0   x1.25   -              refuse:loop start C[1] is not its trial T[2] -- pass 0 runs at the start
//   C  1             INFO   CF>0  x1.25   -              -- the passes before the recorded one are not in the record
//   T  0-1           IN     CF    x1.25   -              -- the trial strain
//   T  2             IN     CF    x1.25   C.1@CF0
//   T  3-5           IN     CF    x1.25   -
//   T  0-1           IN     R     x1.25   -
//   T  2             OUT    R     x1.25   -              fail: T2end( -- the eps_22 the loop left in the trial strain
//   T  3-5           IN     R     x1.25   -
//   O  0-1           OUT    all   +1      -              refuse:disagree with the O line -- status pieces; PIN to K status= subInfo=
//   O  2-16          OUT    all   x1.25   -              -- kappa omega_t omega_c sig[6] sigEff[6]
//   D  0-35          OUT    all   x1.25   -              -- the 6x6 tangent the pyd dumped (C records)
//   M  0-35          OUT    all   x1.25   -              -- the same for F and R (last pass), testbed-recorded
//   W  0-11          OUT    all   x1.25   -              -- eig3sym of the effective stress: the libm canary
//   Z  0             PIN    all   +1      -              refuse:is not its position -- per pass of an R record: the pass index
//   Z  1-3           OUT    all   x1.25   -              -- eps22 s22 D22
//   Z  4-5           OUT    all   +1      -              -- status pieces
//   N  0-39          OUT    all   x1.25   -              -- the full State returnMap returned (note 3)
//   N  40            PIN    all   +1      -              refuse:N line's subInfo is not the K line's -- subInfo; PIN to K subInfo=
//   ---- END FIELD TABLE ----
// Notes to the table:
//   1. K cls= is the dump's return class (0 direct, 1 sub-incremented, 2 floor A, 3 floor B, 4 floor C, 5 direct
//      failure), DERIVED from the replayed returnMap: status 0 -> 0 if subInfo is 0, else 1; status != 0 -> 4 if
//      maxSubIncr > 0, else 5. This tree's kernel has no class API (the probe's DiagMap is not in it): the return status
//      and State::subInfo are its whole report, and they determine every class the probe could write. Class 3 is never
//      assigned by the probe; class 2 needs a piece f < 2^-maxSubIncr, impossible because every piece and every partial
//      sum is an exact multiple of 2^-maxSubIncr; a failure with maxSubIncr > 0 leaves only through floor C. K subInfo=
//      is -1 when status != 0, else 0 or (with maxSubIncr > 0) the piece count >= 2.
//   2. #@ binary= is opensees.pyd@ec969aef2 for the libmmt fixture; `record` keeps it only when it reproduced every
//      O/D/W/Z value and K outcome field on the same libm, otherwise it writes concrete3d_replay-record.
//   3. N: eps[6] sig[6] sigEff[6] kp et_max kdt1 kdt2 kdc kdc1 kdc2 sigtMax sigcMax eqc etPrev wt wc dwt dwc dt_n
//      depl[6] subInfo, 41 values, every State field, what the next step would read (R: after the last pass); required
//      on every record. Testbed-recorded under the fixture's libm, like M. On the libmmt fixture they were recorded by a
//      build that reproduced every pyd-written value of every record bit for bit.
//   4. H: in.etPrev reaches N.eqc and the outputs through the CDPM2 compressive history (compDrive = 1 on every record);
//      the kernel sets N.etPrev from this step's drive. I: wt wc reach N.dwt N.dwc. P eta and implex are mode
//      switches, 0 on every record, so the Duvaut-Lions and IMPL-EX paths are a coverage limit of this fixture.
//   5. K kind= is IN, but every changed kind is refused at load by its line shape (Z lines, the D/M tangent letter):
//      the R rows (C[1]@R, T[i]@R) prove the R replay reads the loop start and strain the C/F replay does not.
//   6. A C/F record's C[1] at pass > 0 is INFO: the passes before the recorded one are not in the record, so the path
//      from C[1] to T[2] is not replayed. `record` writes its banner as #= lines and drops a source's #= lines.
// Tangent lines: a C record must carry D, an R or F record M.
// The sigma_22 loop below is a COPY of the wrapper's: the fixture's R records prove the copy matches the wrapper the
// records were written by; a later change to LadrunoConcrete3D::setTrialStrain must be mirrored here.
#include "SRC/material/nD/LadrunoConcrete3DKernel.h"
#include <algorithm>
#include <cerrno>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
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

// A refusal at load (exit 2). It is thrown, not exited, so `contract` can run a refused mutant in-process; main()
// prints it as "REPLAY ERROR: ..." and exits 2, exactly as `check` always has.
struct Refusal { std::string msg; };
[[noreturn]] static void die(const std::string& m) { throw Refusal{m}; }
static std::string strf(const char* f, ...) {
    va_list a; va_start(a, f); char b[4096]; std::vsnprintf(b, sizeof b, f, a); va_end(a); return b;
}

static const int N_STATE = 41;
static const char* nName[N_STATE] = {
    "eps0", "eps1", "eps2", "eps3", "eps4", "eps5", "sig0", "sig1", "sig2", "sig3", "sig4", "sig5",
    "sigEff0", "sigEff1", "sigEff2", "sigEff3", "sigEff4", "sigEff5", "kp", "et_max", "kdt1", "kdt2", "kdc", "kdc1",
    "kdc2", "sigtMax", "sigcMax", "eqc", "etPrev", "wt", "wc", "dwt", "dwc", "dt_n",
    "depl0", "depl1", "depl2", "depl3", "depl4", "depl5", "subInfo"};
static const char* oName[17] = {"status", "pieces", "kappa", "omega_t", "omega_c", "sig0", "sig1", "sig2", "sig3", "sig4",
                                "sig5", "sigEff0", "sigEff1", "sigEff2", "sigEff3", "sigEff4", "sigEff5"};
static const char* hName[11] = {"kp", "et_max", "kdt1", "kdt2", "kdc", "kdc1", "kdc2", "sigtMax", "sigcMax", "eqc", "etPrev"};
static const char* iName[11] = {"wt", "wc", "dwt", "dwc", "dt_n", "depl0", "depl1", "depl2", "depl3", "depl4", "depl5"};
static const char* zName[6] = {"index", "eps22", "s22", "D22", "status", "pieces"};

// ---- THE FIELD TABLE: the single source of the FIELD CLASSIFICATION -------------------------------------------------
// One entry per field, or per run of columns that share a class and a scope. It drives:
//   * the loader: each numeric line's column count, each key=value line's key set, every enumerated domain;
//   * `fields`: the table printed; the header block between BEGIN/END FIELD TABLE is that output verbatim, and
//     `fields --check <this .cpp>` (run by ci.sh and the .bat) fails when the two differ;
//   * `contract`: one generated row per column (per domain value for {}), its perturbation and its expected outcome.
// So a new column cannot be loaded without an entry, and an entry is a row: a field cannot be added unclassified,
// and a class cannot be wrong without its row going red.
//   cols    a key (key=value lines) or an index range a-b (numeric lines; Q counts its key=value tokens)
//   class   IN | OUT | PIN | UNREAD | INFO (the header's definitions)
//   scope   the records a row edits: all | R | CF | CF0 (C/F at pass 0) | CF>0 | dt>0 | dt<0 | dt<=0 (the record's
//           P dt) | cls0 (K cls=0); lines before the first K line are always `all`
//   pert    x1.25 (a double: v*1.25, 0 -> 0.25) | +1 (an integer) | =V (set V) | {} (one row per domain value, =V; a
//           value no record in scope lacks is reported as skipped) | ~last (flip the last hex digit) | - (no row)
//   tie     edits the same row makes elsewhere, a pin it must keep: <line>.<col>[@scope][=V], comma-separated;
//           without =V the row's own perturbation, applied to that field's own (bit-equal) value
//   domain  the values the loader accepts (comma-separated), "" when not enumerated
//   expect  "" = the class's outcome: IN fails (exit 1, the first mismatch names a compared field); OUT fails naming
//           its own compare; UNREAD and INFO pass. refuse:<text> = exit 2, refused at load with <text> (a pinned
//           interaction); fail:<text> = exit 1 with <text> in the report.
struct FieldSpec { const char *line, *cols, *cls, *scope, *pert, *tie, *domain, *expect, *note; };
static const FieldSpec FIELDS[] = {
    {"#",  "text",          "INFO",   "all",  "-",     "", "", "", "header comments: only before the first K line"},
    {"#@", "libm",          "PIN",    "all",  "=other", "", "", "refuse:do not replay bit-exactly under", "the outputs' libm; --cross waives"},
    {"#@", "kernel-blob",   "OUT",    "all",  "~last", "", "", "fail:KERNEL-BLOB MISMATCH", "the kernel header blob, vs --kernel-blob"},
    {"#@", "binary",        "INFO",   "all",  "=x",    "", "", "", "who wrote the O/D/W/Z values and the K outcome fields (note 2)"},
    {"#@", "M-lines",       "INFO",   "all",  "=x",    "", "", "", "who wrote the M lines"},
    {"#@", "N-lines",       "INFO",   "all",  "=x",    "", "", "", "who wrote the N lines"},
    {"Q",  "0-20",          "OUT",    "all",  "x1.25", "", "", "fail:libm probes (Q line, 21 values): DIFFER", "21 libm probes, vs this build's"},
    {"K",  "kind",          "IN",     "R",    "=F",    "", "C,F,R", "refuse:Z lines on a non-R record", "C/F: one integrate, R: the sigma_22 loop"},
    {"K",  "seq",           "INFO",   "all",  "+1",    "", "", "", "the dump's identifiers and counters: the kernel has none"},
    {"K",  "tag",           "INFO",   "all",  "+1",    "", "", "", "(as seq)"},
    {"K",  "ele",           "INFO",   "all",  "+1",    "", "", "", "(as seq)"},
    {"K",  "calls",         "INFO",   "all",  "+1",    "", "", "", "(as seq)"},
    {"K",  "integ",         "INFO",   "all",  "+1",    "", "", "", "(as seq)"},
    {"K",  "dim",           "INFO",   "CF",   "+1",    "", "", "", "the wrapper's element dimension: the kernel is always 3-D"},
    {"K",  "dim",           "PIN",    "R",    "=5",    "", "", "refuse:without a condensed sigma_22 trace", "a condensed element: 3 or 4"},
    {"K",  "pass",          "OUT",    "R",    "+1",    "", "", "fail: reason ", "the rule-B refusal code of the replayed loop"},
    {"K",  "pass",          "INFO",   "CF",   "+1",    "", "", "", "the sigma_22 pass of the recorded integrate(); pass 0 ties C[1] to T[2]"},
    {"K",  "status",        "PIN",    "all",  "+1",    "", "", "refuse:disagree with the O line", "= O[0]"},
    {"K",  "subInfo",       "PIN",    "all",  "+1",    "", "", "refuse:disagree with the O line", "= O[1] and N[40], in its domain (note 1)"},
    {"K",  "cls",           "OUT",    "all",  "+1",    "", "", "refuse:is not what the dump writes", "derived (note 1); its load pin refuses first"},
    {"G",  "dump",          "INFO",   "all",  "=unset", "", "unset,COMPATIBLE", "", "why the record was picked"},
    {"G",  "class",         "INFO",   "all",  "=direct", "", "refusal,floor-C,sub-increment,kink,direct", "", "groups the report"},
    {"G",  "pick",          "INFO",   "all",  "=x",    "", "", "", "the selection rule"},
    {"P",  "E",             "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "nu",            "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "rho",           "UNREAD", "all",  "x1.25", "", "", "", "the kernel has no use of Params::rho"},
    {"P",  "fc",            "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "ft",            "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "e",             "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "m0",            "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "Gf",            "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "Gc",            "UNREAD", "all",  "x1.25", "", "", "", "read only when epsFc <= 0 (epsFcOf); epsFc > 0 on every record"},
    {"P",  "Df",            "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "As",            "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "qh0",           "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "Hp",            "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "Ah",            "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "Bh",            "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "Ch",            "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "Dh",            "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "ctTemper",      "IN",     "all",  "{}",    "", "0,1,2", "", "-ctTemper none|alphat|proj"},
    {"P",  "tcTemper",      "IN",     "all",  "{}",    "", "0,2", "", "-tcTemper none|proj"},
    {"P",  "tensionLaw",    "IN",     "all",  "{}",    "", "0,1", "", "-tensionLaw exp|bilinear"},
    {"P",  "epsFc",         "IN",     "all",  "x1.25", "", "", "", ""},
    {"P",  "flowPotential", "IN",     "cls0", "=0",    "P.maxSubIncr=0", "0,1", "", "cls 0 only: elsewhere the cls pin refuses it"},
    {"P",  "maxSubIncr",    "PIN",    "all",  "=0",    "", "", "refuse:is not a value the wrapper writes", "10 iff flowPotential=1, else 0"},
    {"P",  "compDrive",     "IN",     "all",  "{}",    "", "0,1", "", "-compressionDrive legacy|cdpm2"},
    {"P",  "eta",           "IN",     "all",  "x1.25", "", "", "", "0 on every record: 0 -> 0.25 switches Duvaut-Lions on"},
    {"P",  "implex",        "IN",     "all",  "{}",    "", "0,1", "", "0 on every record: the IMPL-EX path is a coverage limit"},
    {"P",  "implexRmax",    "PIN",    "all",  "x1.25", "", "", "refuse:is not the dump's constant 2", "the dump's constant; IMPL-EX only"},
    {"P",  "lch",           "IN",     "all",  "x1.25", "P.lch_ref", "", "", ""},
    {"P",  "lch_ref",       "PIN",    "all",  "x1.25", "", "", "refuse:is not lch bit for bit", "the dump prints lch twice; never read"},
    {"P",  "ftol",          "PIN",    "all",  "x1.25", "", "", "refuse:this kernel has another", "= this kernel's F_TOL_HONEST (a constant)"},
    {"P",  "dt",            "IN",     "dt>0", "x1.25", "", "", "", "into N.dt_n"},
    {"P",  "dt",            "UNREAD", "dt<0", "x1.25", "", "", "", "only its sign is read (Duvaut-Lions and IMPL-EX off)"},
    {"E",  "0-1",           "IN",     "all",  "x1.25", "", "", "", "committed strain State::eps"},
    {"E",  "2",             "IN",     "all",  "x1.25", "C.0", "", "", ""},
    {"E",  "3-5",           "IN",     "all",  "x1.25", "", "", "", ""},
    {"S",  "0-5",           "UNREAD", "all",  "x1.25", "", "", "", "State::sig: returnMap never reads in.sig (it returns from X)"},
    {"X",  "0-5",           "IN",     "all",  "x1.25", "", "", "", "committed effective stress State::sigEff"},
    {"H",  "0-10",          "IN",     "all",  "x1.25", "", "", "", "kp et_max kdt1 kdt2 kdc kdc1 kdc2 sigtMax sigcMax eqc etPrev"},
    {"I",  "0-1",           "IN",     "all",  "x1.25", "", "", "", "wt wc, into N.dwt N.dwc"},
    {"I",  "2-3",           "UNREAD", "all",  "x1.25", "", "", "", "dwt dwc: IMPL-EX only"},
    {"I",  "4",             "IN",     "dt<=0", "x1.25", "", "", "", "dt_n, into N.dt_n when dt <= 0"},
    {"I",  "4",             "UNREAD", "dt>0", "x1.25", "", "", "", "dt_n: N.dt_n = dt when dt > 0"},
    {"I",  "5-10",          "UNREAD", "all",  "x1.25", "", "", "", "depl[6]: IMPL-EX only"},
    {"C",  "0",             "PIN",    "all",  "x1.25", "", "", "refuse:is not E[2] bit for bit", "the committed eps_22 = E[2]"},
    {"C",  "1",             "IN",     "R",    "x1.25", "", "", "", "the sigma_22 loop's start"},
    {"C",  "1",             "PIN",    "CF0",  "x1.25", "", "", "refuse:loop start C[1] is not its trial T[2]", "pass 0 runs at the start"},
    {"C",  "1",             "INFO",   "CF>0", "x1.25", "", "", "", "the passes before the recorded one are not in the record"},
    {"T",  "0-1",           "IN",     "CF",   "x1.25", "", "", "", "the trial strain"},
    {"T",  "2",             "IN",     "CF",   "x1.25", "C.1@CF0", "", "", ""},
    {"T",  "3-5",           "IN",     "CF",   "x1.25", "", "", "", ""},
    {"T",  "0-1",           "IN",     "R",    "x1.25", "", "", "", ""},
    {"T",  "2",             "OUT",    "R",    "x1.25", "", "", "fail: T2end(", "the eps_22 the loop left in the trial strain"},
    {"T",  "3-5",           "IN",     "R",    "x1.25", "", "", "", ""},
    {"O",  "0-1",           "OUT",    "all",  "+1",    "", "", "refuse:disagree with the O line", "status pieces; PIN to K status= subInfo="},
    {"O",  "2-16",          "OUT",    "all",  "x1.25", "", "", "", "kappa omega_t omega_c sig[6] sigEff[6]"},
    {"D",  "0-35",          "OUT",    "all",  "x1.25", "", "", "", "the 6x6 tangent the pyd dumped (C records)"},
    {"M",  "0-35",          "OUT",    "all",  "x1.25", "", "", "", "the same for F and R (last pass), testbed-recorded"},
    {"W",  "0-11",          "OUT",    "all",  "x1.25", "", "", "", "eig3sym of the effective stress: the libm canary"},
    {"Z",  "0",             "PIN",    "all",  "+1",    "", "", "refuse:is not its position", "per pass of an R record: the pass index"},
    {"Z",  "1-3",           "OUT",    "all",  "x1.25", "", "", "", "eps22 s22 D22"},
    {"Z",  "4-5",           "OUT",    "all",  "+1",    "", "", "", "status pieces"},
    {"N",  "0-39",          "OUT",    "all",  "x1.25", "", "", "", "the full State returnMap returned (note 3)"},
    {"N",  "40",            "PIN",    "all",  "+1",    "", "", "refuse:N line's subInfo is not the K line's", "subInfo; PIN to K subInfo="},
};
static const char* LINE_KINDS[] = {"#", "#@", "Q", "K", "G", "P", "E", "S", "X", "H", "I", "C", "T", "O", "D", "M", "W", "Z", "N"};

static bool indexed(const FieldSpec& f) { return f.cols[0] >= '0' && f.cols[0] <= '9'; }
static void colRange(const FieldSpec& f, int& a, int& b) {
    const char* dash = std::strchr(f.cols, '-'); a = std::atoi(f.cols); b = dash ? std::atoi(dash + 1) : a;
}
static int ncols(const char* line) {   // a numeric line's column count, from the table
    int n = 0; for (const FieldSpec& f : FIELDS) if (!std::strcmp(f.line, line) && indexed(f)) { int a, b; colRange(f, a, b); n = std::max(n, b + 1); }
    return n;
}
static std::set<std::string> keysOf(const char* line) {
    std::set<std::string> s; for (const FieldSpec& f : FIELDS) if (!std::strcmp(f.line, line) && !indexed(f)) s.insert(f.cols);
    return s;
}
static std::vector<std::string> split(const std::string& s, char c) {
    std::vector<std::string> v; if (s.empty()) return v; size_t p = 0, q;
    while ((q = s.find(c, p)) != std::string::npos) { v.push_back(s.substr(p, q - p)); p = q + 1; }
    v.push_back(s.substr(p)); return v;
}
static std::vector<std::string> domainOf(const char* line, const char* key) {
    for (const FieldSpec& f : FIELDS) if (!std::strcmp(f.line, line) && !std::strcmp(f.cols, key) && *f.domain) return split(f.domain, ',');
    die(std::string("field table: ") + line + "." + key + " has no domain");
}
// the label of numeric column i (row ids) and, for a compared line, the name compare() reports it by
static std::string colLabel(const std::string& line, int i) {
    if (line == "H") return hName[i];
    if (line == "I") return iName[i];
    if (line == "O") return oName[i];
    if (line == "N") return nName[i];
    if (line == "Z") return zName[i];
    return "";
}
static std::string colName(const std::string& line, int i) {
    if (line == "O") return oName[i];
    if (line == "N") return std::string("N.") + nName[i];
    if (line == "D" || line == "M") return line + std::to_string(i / 6) + std::to_string(i % 6);
    if (line == "W") return "eig3sym" + std::to_string(i);
    if (line == "Z") return std::string("pass0.") + zName[i];   // every R record has a pass 0
    return "";
}

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
// Column counts, key sets and domains come from the field table.
static Fixture loadStream(std::istream& in, bool needN = true) {
    Fixture fx; std::string line; Rec* r = 0; bool haveMeta = false;
    static const std::set<std::string> kKeys = keysOf("K"), gKeys = keysOf("G"), metaKeys = keysOf("#@");
    static const std::vector<std::string> kinds = domainOf("K", "kind"), dumps = domainOf("G", "dump"), classes = domainOf("G", "class");
    auto inDomain = [](const std::vector<std::string>& d, const std::string& v) { return std::find(d.begin(), d.end(), v) != d.end(); };
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
            for (const std::string& key : kKeys) if (!r->k.count(key)) die("K line lacks " + key + ": " + line.substr(0, 60));
            for (const auto& e : r->k) if (!kKeys.count(e.first)) die("K line has an unexpected key " + e.first + ": " + line.substr(0, 60));
            const std::string kind = r->k.at("kind");
            if (!inDomain(kinds, kind)) die("unknown record kind '" + kind + "' (" + line.substr(0, 40) + ")");
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
            bool keysOk = r->g.size() == gKeys.size();
            for (const std::string& key : gKeys) keysOk = keysOk && r->g.count(key);
            if (!keysOk) die("G line must carry exactly dump= class= pick=: " + line);
            if (!inDomain(dumps, r->g.at("dump"))) die("G dump= is not unset|COMPATIBLE: " + line);
            if (!inDomain(classes, r->g.at("class"))) die("G class= is not a record class: " + line);
            break; }
        case 'P': {
            r->pLine = line; kv(line, r->pkv);
            static const std::set<std::string> pKeys = keysOf("P");
            for (const std::string& key : pKeys) if (!r->pkv.count(key)) die("P line lacks " + key + ": " + r->kLine);
            for (const auto& e : r->pkv) if (!pKeys.count(e.first)) die("P line has a key this testbed does not map: " + e.first + " (" + r->kLine + ")");
            break; }
        case 'E': copy(nums(line, ncols("E")), r->eps); r->stateLines.push_back(line); break;
        case 'S': copy(nums(line, ncols("S")), r->sig); r->stateLines.push_back(line); break;
        case 'X': copy(nums(line, ncols("X")), r->sigEff); r->stateLines.push_back(line); break;
        case 'H': copy(nums(line, ncols("H")), r->h); r->stateLines.push_back(line); break;
        case 'I': copy(nums(line, ncols("I")), r->im); r->stateLines.push_back(line); break;
        case 'C': copy(nums(line, ncols("C")), r->c); r->stateLines.push_back(line); break;
        case 'T': copy(nums(line, ncols("T")), r->strain); r->tTok = toks(line); r->stateLines.push_back(line); break;
        case 'O': copy(nums(line, ncols("O"), {0, 1}), r->o); r->hasO = true; break;
        case 'D': case 'M': copy(nums(line, ncols(c == 'D' ? "D" : "M")), r->D); r->hasD = true; r->dKind = c; break;
        case 'W': copy(nums(line, ncols("W")), r->w); r->hasW = true; break;
        case 'N': copy(nums(line, ncols("N"), {N_STATE - 1}), r->n); r->hasN = true; break;
        case 'Z': { auto v = nums(line, ncols("Z"), {0, 4, 5}); Pass p{(int)v[0], (int)v[4], (int)v[5], v[1], v[2], v[3]};
                    if (p.p != (int)r->z.size()) die("Z pass index " + std::to_string(p.p) + " is not its position " + std::to_string(r->z.size()) + ": " + r->kLine);
                    r->z.push_back(p); break; }
        default: die(std::string("unknown line kind '") + c + "' (the fixture is the record-dump subset: Q K G P E S X H I C T O D M W Z N)");
        }
    }
    if (fx.recs.empty()) die("the fixture has no records");
    if (fx.qLine.empty()) die("the fixture has no Q line (libm probes)");
    for (const char* key : {"libm", "kernel-blob"}) if (!fx.meta.count(key)) die(std::string("the #@ line lacks ") + key);
    for (const auto& e : fx.meta) if (!metaKeys.count(e.first)) die("the #@ line has an unexpected key " + e.first);
    { const std::string& b = fx.meta.at("kernel-blob");
      if (b.size() != 40 || b.find_first_not_of("0123456789abcdef") != std::string::npos) die("the #@ kernel-blob is not a 40-hex git blob: " + b); }
    for (const Rec& q : fx.recs) {
        for (char c : std::string(needN ? "GPESXHICTOWN" : "GPESXHICTOW")) if (!q.seen.count(c)) die(std::string("record lacks its ") + c + " line: " + q.kLine);
        // the line shape first (note 5): a record whose kind= disagrees with its Z lines is refused for that, before any
        // value pin below can refuse it for a pass-dependent reason (a re-recorded R record's rule-B code can be 0)
        if (q.kind != 'R' && !q.z.empty()) die("Z lines on a non-R record: " + q.kLine);
        if (std::memcmp(&q.c[0], &q.eps[2], sizeof(double)) != 0) die("C line's committed eps_22 is not E[2] bit for bit: " + q.kLine);
        if (q.kind != 'R' && q.code == 0 && std::memcmp(&q.c[1], &q.strain[2], sizeof(double)) != 0)
            die("a pass-0 C/F record's loop start C[1] is not its trial T[2] bit for bit: " + q.kLine);
        K::Params pp; double dt = 0; fillParams(q, pp, dt);   // every P line is validated before anything replays
        if (!q.hasD || q.dKind != (q.kind == 'C' ? 'D' : 'M'))
            die(std::string("a ") + q.kind + " record must carry its 6x6 tangent as a " + (q.kind == 'C' ? "D" : "M") + " line: " + q.kLine);
        if (q.kind == 'R' && (q.z.empty() || (q.k.at("dim") != "4" && q.k.at("dim") != "3"))) die("R record without a condensed sigma_22 trace: " + q.kLine);
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
static Fixture load(const char* path, bool needN = true) {
    std::ifstream in(path); if (!in) die(std::string("cannot open ") + path);
    return loadStream(in, needN);
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
    // only values LadrunoConcrete3D's option parser can produce (the field table's domains, and maxSubIncr = 10 iff
    // the CDPM2 potential): anything else would replay as some other mode without a word
    auto bad = [&](const char* key, int v) { die(std::string("P ") + key + "=" + std::to_string(v) + " is not a value the wrapper writes: " + r.kLine); };
    auto in = [&](const char* key, int v) {
        for (const std::string& o : domainOf("P", key)) if (std::to_string(v) == o) return;
        bad(key, v);
    };
    in("ctTemper", p.ctTemper); in("tcTemper", p.tcTemper); in("tensionLaw", p.tensionLaw);
    in("flowPotential", p.flowPotential); in("compDrive", p.compDrive); in("implex", implex);
    if (p.maxSubIncr != (p.flowPotential == 1 ? 10 : 0)) bad("maxSubIncr", p.maxSubIncr);
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
    for (int i = 0; i < 17; ++i) { if (i < 2) cmpInt(colName("O", i), (long)x.o[i], (long)r.o[i]); else cmp(colName("O", i), x.o[i], r.o[i]); }
    cmpInt("cls", x.cls, r.cls);
    for (int i = 0; i < 36; ++i) cmp(colName(std::string(1, r.dKind), i), (&x.ev.D[0][0])[i], r.D[i]);
    for (int i = 0; i < 12; ++i) cmp(colName("W", i), x.w[i], r.w[i]);
    if (withN) for (int i = 0; i < N_STATE; ++i) {
        if (i == N_STATE - 1) cmpInt(colName("N", i), (long)x.n[i], (long)r.n[i]);
        else cmp(colName("N", i), x.n[i], r.n[i]);
    }
    if (r.kind == 'R') {
        cmpInt("reason", x.loop.code, r.code);
        cmpInt("passes", (long)x.loop.z.size(), (long)r.z.size());
        cmp("T2end", x.loop.eps22End, r.strain[2]);
        for (size_t p = 0; p < x.loop.z.size() && p < r.z.size(); ++p) {
            const Pass &a = x.loop.z[p], &b = r.z[p]; const std::string at = "pass" + std::to_string(p) + ".";
            cmpInt(at + zName[4], a.status, b.status); cmpInt(at + zName[5], a.sub, b.sub);
            cmp(at + zName[1], a.eps22, b.eps22); cmp(at + zName[2], a.s22, b.s22); cmp(at + zName[3], a.d22, b.d22);
        }
    }
    return nDiff;
}

static void checkBlobArg(const std::string& blob) {
    if (blob.size() != 40 || blob.find_first_not_of("0123456789abcdef") != std::string::npos)
        die("--kernel-blob must be the 40-hex git blob of the kernel header this build compiled (git hash-object)");
}

// The replay is a pure function of a record's inputs: its kind, its P line and its E S X H I C T lines (verbatim, the
// unread ones included). `contract` replays the base fixture once and reuses a replay for a record whose inputs a row
// did not touch; a record it touched is replayed afresh.
static std::string inputKey(const Rec& r) {
    std::string k(1, r.kind); k += "\n" + r.pLine;
    for (const std::string& s : r.stateLines) k += "\n" + s;
    return k;
}
struct ReplayCache { std::map<std::string, Replay> byKey; std::vector<size_t> order; };   // order: cheapest record first

// `check`'s whole verdict, written through `emit`. cache: reuse the base replays; stopAtFirst: stop at the first
// mismatching record (cheapest first), for a row that must fail. Returns 0 pass, 1 fail; a refusal throws.
static int runCheck(const Fixture& fx, bool cross, const std::string& blob, const std::function<void(const std::string&)>& emit,
                    const ReplayCache* cache = 0, bool stopAtFirst = false) {
    const std::string& libm = fx.meta.at("libm");
    emit(strf("fixture: %zu records, recorded with libm=%s kernel-blob=%s; this build: libm=%s kernel-blob=%s\n", fx.recs.size(),
              libm.c_str(), fx.meta.at("kernel-blob").c_str(), C3D_REPLAY_LIBM, blob.c_str()));
    if (libm != C3D_REPLAY_LIBM && !cross)
        die("records written under libm=" + libm + " do not replay bit-exactly under " + C3D_REPLAY_LIBM +
            " (run `record` for record-and-compare on this libm, or --cross to measure the difference)");
    if (cross) emit("CROSS-LIBM comparison (diagnostic): a mismatch here is the libm, not the kernel\n");
    const bool blobOk = (fx.meta.at("kernel-blob") == blob);
    if (!blobOk) emit(strf("KERNEL-BLOB MISMATCH: the fixture was recorded on kernel-blob %s, this build compiled %s: re-record (`record`)\n",
                           fx.meta.at("kernel-blob").c_str(), blob.c_str()));
    int bad = 0, seen = 0; const bool qOk = (qProbe() == fx.qLine);
    emit(strf("libm probes (Q line, 21 values): %s\n", qOk ? "bit-identical" : "DIFFER"));
    std::map<std::string, std::pair<int, int>> byGroup;   // selection class -> (records, bit-exact)
    std::vector<size_t> order;
    if (cache && cache->order.size() == fx.recs.size()) order = cache->order;
    else for (size_t i = 0; i < fx.recs.size(); ++i) order.push_back(i);
    for (size_t ri : order) {
        if (stopAtFirst && (bad || !qOk || !blobOk)) break;
        const Rec& r = fx.recs[ri]; ++seen;
        std::map<std::string, Replay>::const_iterator hit;
        const bool cached = cache && (hit = cache->byKey.find(inputKey(r))) != cache->byKey.end();
        std::string why; const int nDiff = compare(r, cached ? hit->second : replay(r), why);
        const std::map<std::string, std::string>& g = r.g;
        auto& cnt = byGroup[g.at("class")]; ++cnt.first;
        if (!nDiff) ++cnt.second;
        else { ++bad; emit(strf("  MISMATCH %s seq=%s (%s %s): %d values differ:%s%s\n", r.k.at("kind").c_str(), r.k.at("seq").c_str(),
                                g.at("dump").c_str(), g.at("class").c_str(), nDiff, why.c_str(), nDiff > 6 ? " ..." : "")); }
    }
    for (const auto& kvp : byGroup) emit(strf("  %-14s %2d records, bit-exact %2d\n", kvp.first.c_str(), kvp.second.first, kvp.second.second));
    const bool fail = bad || !qOk || !blobOk;
    emit(strf("REPLAY %s: %d of %zu records bit-exact%s%s%s\n", fail ? "FAIL" : "PASS", seen - bad, fx.recs.size(),
              qOk ? "" : ", libm probes differ", blobOk ? "" : ", kernel-blob differs",
              seen < (int)fx.recs.size() ? " (stopped at the first mismatch)" : ""));
    return fail ? 1 : 0;
}
static int modeCheck(const Fixture& fx, bool cross, const std::string& blob) {
    return runCheck(fx, cross, blob, [](const std::string& s) { std::fputs(s.c_str(), stdout); });
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

// ---- fields: the table, printed; checked against the header block -------------------------------------------------
static std::vector<std::string> fieldsText() {
    std::vector<std::string> v;
    v.push_back(strf("%-2s %-13s %-6s %-5s %-7s %-14s %s", "ln", "cols", "class", "scope", "pert", "tie", "expect -- note"));
    for (const FieldSpec& f : FIELDS) {
        std::string tail = *f.expect ? f.expect : "";
        std::string note = f.note; if (*f.domain) note += std::string(note.empty() ? "" : "; ") + "domain " + f.domain;
        if (!note.empty()) tail += (tail.empty() ? "-- " : " -- ") + note;
        std::string s = strf("%-2s %-13s %-6s %-5s %-7s %-14s %s", f.line, f.cols, f.cls, f.scope, f.pert, *f.tie ? f.tie : "-", tail.c_str());
        while (!s.empty() && s.back() == ' ') s.pop_back();
        v.push_back(s);
    }
    return v;
}
static const char* BEGIN_MARK = "BEGIN FIELD TABLE";
static const char* END_MARK = "END FIELD TABLE";
// The block between the markers in <source>, each line "//   " + a line of fieldsText(), must equal fieldsText().
static int modeFieldsCheck(const char* src) {
    std::ifstream in(src); if (!in) die(std::string("cannot open ") + src);
    std::vector<std::string> got; std::string line; int state = 0;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.rfind(std::string("//   ---- ") + BEGIN_MARK, 0) == 0) { if (state) die("two BEGIN FIELD TABLE markers"); state = 1; continue; }
        if (line.rfind(std::string("//   ---- ") + END_MARK, 0) == 0) { if (state != 1) die("END FIELD TABLE without BEGIN"); state = 2; continue; }
        if (state == 1) {
            if (line.rfind("//   ", 0) != 0) die("a field-table header line does not start with '//   ': " + line.substr(0, 60));
            got.push_back(line.substr(5));
        }
    }
    if (state != 2) die(std::string("no BEGIN/END FIELD TABLE block in ") + src);
    const std::vector<std::string> want = fieldsText();
    for (size_t i = 0; i < std::max(got.size(), want.size()); ++i) {
        const std::string a = i < got.size() ? got[i] : "(missing)", b = i < want.size() ? want[i] : "(missing)";
        if (a != b) {
            std::printf("FIELD TABLE DIFFERS from the compiled table at header line %zu:\n  header: %s\n  table:  %s\n"
                        "regenerate the block with `concrete3d_replay fields --header`\n", i + 1, a.c_str(), b.c_str());
            return 1;
        }
    }
    std::printf("FIELD TABLE MATCHES: the header block is the compiled table (%zu entries)\n", want.size() - 1);
    return 0;
}

// ---- the table's own consistency, checked at start-up (every mode) --------------------------------------------------
static const std::set<std::string> SCOPES = {"all", "R", "CF", "CF0", "CF>0", "dt>0", "dt<0", "dt<=0", "cls0"};
struct Tie { std::string line, col, scope, set; };
static std::vector<Tie> parseTies(const FieldSpec& f) {
    std::vector<Tie> v;
    for (const std::string& t : split(f.tie, ',')) {
        Tie x; std::string s = t; size_t p;
        if ((p = s.find('=')) != std::string::npos) { x.set = s.substr(p + 1); s = s.substr(0, p); }
        if ((p = s.find('@')) != std::string::npos) { x.scope = s.substr(p + 1); s = s.substr(0, p); } else x.scope = f.scope;
        if ((p = s.find('.')) == std::string::npos) die(std::string("field table: tie '") + t + "' is not <line>.<col>");
        x.line = s.substr(0, p); x.col = s.substr(p + 1);
        if (!SCOPES.count(x.scope)) die("field table: tie '" + t + "' has an unknown scope");
        v.push_back(x);
    }
    return v;
}
static void tableSelfCheck() {
    std::set<std::string> lines(std::begin(LINE_KINDS), std::end(LINE_KINDS)), used, entries;
    for (const FieldSpec& f : FIELDS) {
        const std::string id = std::string(f.line) + "." + f.cols + "@" + f.scope, c = f.cls, p = f.pert, e = f.expect;
        if (!lines.count(f.line)) die("field table: unknown line kind " + id);
        if (!entries.insert(id).second) die("field table: two entries for " + id);
        used.insert(f.line);
        if (c != "IN" && c != "OUT" && c != "PIN" && c != "UNREAD" && c != "INFO") die("field table: unknown class " + c + " (" + id + ")");
        if (!SCOPES.count(f.scope)) die("field table: unknown scope (" + id + ")");
        if (p != "-" && p != "x1.25" && p != "+1" && p != "~last" && p != "{}" && p[0] != '=') die("field table: unknown pert " + p + " (" + id + ")");
        if (p == "{}" && !*f.domain) die("field table: pert {} needs a domain (" + id + ")");
        if (c == "PIN" && e.rfind("refuse:", 0) != 0) die("field table: a PIN entry must name its refusal (" + id + ")");
        if ((c == "UNREAD" || c == "INFO") && !e.empty()) die("field table: an UNREAD or INFO entry passes; it takes no expect (" + id + ")");
        if (!e.empty() && e.rfind("refuse:", 0) != 0 && e.rfind("fail:", 0) != 0) die("field table: expect is refuse:|fail: (" + id + ")");
        if (c == "OUT" && e.empty() && (!indexed(f) || colName(f.line, 0).empty())) die("field table: an OUT entry without a compare name needs fail: (" + id + ")");
        if (indexed(f) && std::string(f.line) != "Q") {   // only lines whose columns are fed or compared in the code
            int a, b; colRange(f, a, b);
            if (a > b) die("field table: bad range (" + id + ")");
        }
        if (*f.tie) for (const Tie& t : parseTies(f)) if (!lines.count(t.line)) die("field table: tie on an unknown line (" + id + ")");
    }
    for (const std::string& l : lines) if (!used.count(l)) die("field table: no entry for line kind " + l);
    // every column of every numeric line is classified, and the table's counts are the code's
    const std::map<std::string, int> code = {{"E", 6}, {"S", 6}, {"X", 6}, {"H", 11}, {"I", 11}, {"C", 2}, {"T", 6}, {"O", 17},
                                             {"D", 36}, {"M", 36}, {"W", 12}, {"Z", 6}, {"N", N_STATE}, {"Q", (int)toks(qProbe()).size()}};
    for (const auto& e : code) {
        const int n = ncols(e.first.c_str());
        if (n != e.second) die(strf("field table: line %s has %d columns, the code %d", e.first.c_str(), n, e.second));
        for (int i = 0; i < n; ++i) {
            bool cov = false;
            for (const FieldSpec& f : FIELDS) if (e.first == f.line && indexed(f)) { int a, b; colRange(f, a, b); cov = cov || (a <= i && i <= b); }
            if (!cov) die(strf("field table: column %s[%d] is not classified", e.first.c_str(), i));
        }
    }
}

// ---- contract: the rows generated from the table, run in-process --------------------------------------------------
struct Row { std::string id; const FieldSpec* f; std::string col; std::string pert; bool skipNoop; };
static std::vector<Row> rowsOf() {
    std::vector<Row> v;
    for (const FieldSpec& f : FIELDS) {
        if (!std::strcmp(f.pert, "-")) continue;
        std::vector<std::string> perts;
        if (!std::strcmp(f.pert, "{}")) for (const std::string& d : split(f.domain, ',')) perts.push_back("=" + d);
        else perts.push_back(f.pert);
        std::vector<std::string> cols;
        if (indexed(f)) { int a, b; colRange(f, a, b); for (int i = a; i <= b; ++i) cols.push_back(std::to_string(i)); }
        else cols.push_back(f.cols);
        for (const std::string& c : cols) for (const std::string& p : perts) {
            std::string id = std::string(f.line) + (indexed(f) ? "[" + c + "]" : "." + c);
            if (indexed(f)) { const std::string lab = colLabel(f.line, std::atoi(c.c_str())); if (!lab.empty()) id += "(" + lab + ")"; }
            if (std::strcmp(f.scope, "all")) id += std::string("@") + f.scope;
            if (!std::strcmp(f.pert, "{}")) id += p;
            v.push_back(Row{id, &f, c, p, !std::strcmp(f.pert, "{}")});
        }
    }
    return v;
}
static std::string expectText(const Row& r) {
    const std::string c = r.f->cls, e = r.f->expect;
    if (e.rfind("refuse:", 0) == 0) return "exit 2, refused at load: '" + e.substr(7) + "'";
    if (e.rfind("fail:", 0) == 0) return "exit 1: '" + e.substr(5) + "'";
    if (c == "UNREAD" || c == "INFO") return "exit 0, REPLAY PASS";
    if (c == "OUT") return "exit 1 naming " + colName(r.f->line, std::atoi(r.col.c_str()));
    return "exit 1 naming a compared field";
}
static int modeFieldsRows() {
    for (const Row& r : rowsOf())
        std::printf("%-26s %-6s %-5s %-9s %s\n", r.id.c_str(), r.f->cls, r.f->scope, r.pert.c_str(), expectText(r).c_str());
    return 0;
}

// text-level mutation: the fixture's lines with one row's edit (and its ties), by token, never through the loader
struct RecCtx { std::string kind; int pass = 0, cls = 0; double dt = 0; };
static bool inScope(const std::string& scope, int ri, const std::vector<RecCtx>& ctx) {
    if (scope == "all") return true;
    if (ri < 0) return false;
    const RecCtx& c = ctx[ri];
    if (scope == "R") return c.kind == "R";
    if (scope == "CF") return c.kind != "R";
    if (scope == "CF0") return c.kind != "R" && c.pass == 0;
    if (scope == "CF>0") return c.kind != "R" && c.pass > 0;
    if (scope == "dt>0") return c.dt > 0;
    if (scope == "dt<0") return c.dt < 0;
    if (scope == "dt<=0") return c.dt <= 0;
    if (scope == "cls0") return c.cls == 0;
    die("unknown scope " + scope);
}
static std::string perturb(const std::string& v, const std::string& pert) {
    if (pert == "x1.25") { const double x = std::strtod(v.c_str(), 0); return strf("%.17g", x == 0 ? 0.25 : x * 1.25); }
    if (pert == "+1") return std::to_string(std::strtoll(v.c_str(), 0, 10) + 1);
    if (pert == "~last") { std::string s = v; if (!s.empty()) s.back() = s.back() == '0' ? '1' : '0'; return s; }
    if (!pert.empty() && pert[0] == '=') return pert.substr(1);
    die("cannot apply pert " + pert);
}
static std::string lineKind(const std::string& l) {
    if (l.rfind("#@ ", 0) == 0) return "#@";
    if (!l.empty() && l[0] == '#') return "#";
    return l.substr(0, 1);
}
static int editLines(std::vector<std::string>& L, const std::vector<int>& recOf, const std::vector<RecCtx>& ctx,
                     const std::string& line, const std::string& col, const std::string& scope, const std::string& pert) {
    const bool idx = col[0] >= '0' && col[0] <= '9'; int n = 0;
    for (size_t i = 0; i < L.size(); ++i) {
        if (lineKind(L[i]) != line || !inScope(scope, recOf[i], ctx)) continue;
        std::vector<std::string> t = split(L[i], ' ');
        std::string* tok = 0; size_t eq = 0;
        if (idx) { const size_t k = 1 + std::atoi(col.c_str()); if (k < t.size()) tok = &t[k]; }
        else for (size_t k = 1; k < t.size(); ++k) if (t[k].rfind(col + "=", 0) == 0) { tok = &t[k]; break; }
        if (!tok) die("contract: " + line + " " + col + " is not on the line: " + L[i].substr(0, 60));
        if (!idx || line == "Q") { eq = tok->find('='); if (eq == std::string::npos) die("contract: not key=value: " + *tok); }
        if (!idx || line == "Q") *tok = tok->substr(0, eq + 1) + perturb(tok->substr(eq + 1), pert);
        else *tok = perturb(*tok, pert);
        std::string s = t[0]; for (size_t k = 1; k < t.size(); ++k) s += " " + t[k];
        if (s != L[i]) ++n;
        L[i] = s;
    }
    return n;
}
static std::vector<std::string> readLines(const char* path) {
    std::ifstream in(path, std::ios::binary); if (!in) die(std::string("cannot open ") + path);
    std::vector<std::string> v; std::string l;
    while (std::getline(in, l)) { if (!l.empty() && l.back() == '\r') l.pop_back(); v.push_back(l); }
    return v;
}
static std::string joinLines(const std::vector<std::string>& L) { std::string s; for (const std::string& l : L) s += l + "\n"; return s; }

// The first compared field a MISMATCH line names ("values differ: <name>(..." or "<name> a!=b")
static std::string firstNamed(const std::string& text) {
    const size_t p = text.find("values differ: "); if (p == std::string::npos) return "";
    const size_t a = p + 15, b = text.find_first_of("( \n", a);
    return text.substr(a, b == std::string::npos ? std::string::npos : b - a);
}
static bool isCompareName(const std::string& n) {
    static std::set<std::string> names;
    if (names.empty()) {
        for (const char* l : {"O", "D", "M", "W", "N"}) for (int i = 0; i < ncols(l); ++i) names.insert(colName(l, i));
        for (const char* s : {"cls", "reason", "passes", "T2end"}) names.insert(s);
    }
    if (names.count(n)) return true;
    if (n.rfind("pass", 0) == 0) {   // pass<p>.<Z field>
        const size_t d = n.find('.'); if (d == std::string::npos || d == 4 || n.find_first_not_of("0123456789", 4) != d) return false;
        for (int i = 1; i < 6; ++i) if (n.substr(d + 1) == zName[i]) return true;
    }
    return false;
}

static int modeContract(const char* path, const std::string& blob) {
    const auto t0 = std::chrono::steady_clock::now();
    const std::vector<std::string> base = readLines(path);
    // the base fixture must pass first (the gate); its replays become the cache, cheapest record first
    const Fixture fx = load(path);
    ReplayCache cache; std::vector<std::pair<double, size_t>> cost;
    for (size_t i = 0; i < fx.recs.size(); ++i) {
        const auto a = std::chrono::steady_clock::now();
        cache.byKey[inputKey(fx.recs[i])] = replay(fx.recs[i]);
        cost.push_back({std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count(), i});
    }
    std::sort(cost.begin(), cost.end());
    for (const auto& c : cost) cache.order.push_back(c.second);
    std::string baseText;
    if (runCheck(fx, false, blob, [&](const std::string& s) { baseText += s; }, &cache) != 0) {
        std::fputs(baseText.c_str(), stdout);
        std::printf("CONTRACT ERROR: the base fixture does not replay (the gate fails): no row can be judged\n");
        return 2;
    }
    // per-line record index and each record's scope context, read from the text
    std::vector<int> recOf(base.size(), -1); std::vector<RecCtx> ctx; int ri = -1;
    for (size_t i = 0; i < base.size(); ++i) {
        if (base[i].rfind("K ", 0) == 0) {
            ++ri; ctx.emplace_back(); std::map<std::string, std::string> k; kv(base[i], k);
            ctx.back().kind = k.at("kind"); ctx.back().pass = std::atoi(k.at("pass").c_str()); ctx.back().cls = std::atoi(k.at("cls").c_str());
        } else if (base[i].rfind("P ", 0) == 0 && ri >= 0) {
            std::map<std::string, std::string> p; kv(base[i], p); ctx[ri].dt = std::strtod(p.at("dt").c_str(), 0);
        }
        recOf[i] = ri;
    }
    const std::string nRec = std::to_string(fx.recs.size());
    std::map<std::string, std::pair<int, int>> byCls;   // class -> (rows run, rows as designed)
    int bad = 0, skipped = 0;
    std::printf("contract: %zu rows generated from the field table, on %s (%s records)\n", rowsOf().size(), path, nRec.c_str());
    for (const Row& row : rowsOf()) {
        const FieldSpec& f = *row.f; const std::string cls = f.cls, e = f.expect;
        std::vector<std::string> L = base;
        int n = editLines(L, recOf, ctx, f.line, row.col, f.scope, row.pert);
        if (*f.tie) for (const Tie& t : parseTies(f)) editLines(L, recOf, ctx, t.line, t.col, t.scope, t.set.empty() ? row.pert : "=" + t.set);
        if (n == 0) {
            if (row.skipNoop) { ++skipped; std::printf("  skip %-6s %-26s every record in scope already has %s\n", cls.c_str(), row.id.c_str(), row.pert.c_str() + 1); continue; }
            ++bad; std::printf("CONTRACT FAIL: %s (%s): the row changed nothing (scope %s matches no %s line)\n", row.id.c_str(), cls.c_str(), f.scope, f.line);
            continue;
        }
        const bool mustFail = !(cls == "UNREAD" || cls == "INFO");
        int rc; std::string text;
        try {
            std::istringstream is(joinLines(L));
            const Fixture m = loadStream(is);
            rc = runCheck(m, false, blob, [&](const std::string& s) { text += s; }, &cache, mustFail);
        } catch (const Refusal& x) { rc = 2; text += "REPLAY ERROR: " + x.msg + "\n"; }
        // what the row must give
        int wantRc = mustFail ? 1 : 0; std::string want, got;
        bool ok;
        if (e.rfind("refuse:", 0) == 0) { wantRc = 2; want = e.substr(7); ok = rc == 2 && text.find(want) != std::string::npos; }
        else if (e.rfind("fail:", 0) == 0) { want = e.substr(5); ok = rc == 1 && text.find(want) != std::string::npos; }
        else if (!mustFail) { want = "REPLAY PASS: " + nRec + " of " + nRec; ok = rc == 0 && text.find(want) != std::string::npos; }
        else if (cls == "OUT") {
            want = colName(f.line, std::atoi(row.col.c_str()));
            ok = rc == 1 && (text.find(" " + want + "(") != std::string::npos || text.find(" " + want + " ") != std::string::npos);
        } else { want = "a compared field"; got = firstNamed(text); ok = rc == 1 && isCompareName(got); }
        std::string shown;
        for (const std::string& l : split(text, '\n'))
            if (l.find("MISMATCH") != std::string::npos || l.find("REPLAY ERROR") != std::string::npos || l.find("REPLAY PASS") != std::string::npos) { shown = l; break; }
        if (shown.empty()) for (const std::string& l : split(text, '\n')) if (l.find("REPLAY FAIL") != std::string::npos) shown = l;
        while (!shown.empty() && shown[0] == ' ') shown.erase(0, 1);
        auto& c = byCls[cls]; ++c.first;
        if (ok) { ++c.second; std::printf("  ok   %-6s %-26s exit %d: %s\n", cls.c_str(), row.id.c_str(), rc, shown.substr(0, 110).c_str()); }
        else {
            ++bad;
            std::printf("CONTRACT FAIL: %s (%s, scope %s, %s): expected %s, got exit %d: %s\n", row.id.c_str(), cls.c_str(), f.scope,
                        row.pert.c_str(), (wantRc == 2 ? "exit 2 with '" + want + "'" : "exit " + std::to_string(wantRc) + " " + want).c_str(),
                        rc, shown.substr(0, 160).c_str());
        }
    }
    std::string sum;
    for (const auto& c : byCls) sum += strf(" %s %d/%d", c.first.c_str(), c.second.second, c.second.first);
    std::printf("CONTRACT %s: rows as designed by class:%s; %d skipped (domain value already everywhere); %.1f s\n",
                bad ? "FAIL" : "PASS", sum.c_str(), skipped, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
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

static int run(int argc, char** argv) {
    checkRuntimeLibm();
    tableSelfCheck();
    if (argc == 2 && std::strcmp(argv[1], "fields") == 0) { for (const std::string& s : fieldsText()) std::printf("%s\n", s.c_str()); return 0; }
    if (argc == 3 && std::strcmp(argv[1], "fields") == 0 && std::strcmp(argv[2], "--header") == 0) {
        std::printf("//   ---- %s (the output of `concrete3d_replay fields`; `fields --check` compares) ----\n", BEGIN_MARK);
        for (const std::string& s : fieldsText()) std::printf("//   %s\n", s.c_str());
        std::printf("//   ---- %s ----\n", END_MARK);
        return 0;
    }
    if (argc == 3 && std::strcmp(argv[1], "fields") == 0 && std::strcmp(argv[2], "--rows") == 0) return modeFieldsRows();
    if (argc == 4 && std::strcmp(argv[1], "fields") == 0 && std::strcmp(argv[2], "--check") == 0) return modeFieldsCheck(argv[3]);
    if (argc == 5 && std::strcmp(argv[1], "contract") == 0 && std::strcmp(argv[3], "--kernel-blob") == 0) {
        checkBlobArg(argv[4]);
        return modeContract(argv[2], argv[4]);
    }
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
    std::fprintf(stderr, "usage: %s check <fixture> --kernel-blob <sha> [--cross] | record <fixture> <out> --kernel-blob <sha> | "
                         "inputs-match <a> <b> | fields [--header|--rows|--check <concrete3d_replay.cpp>] | contract <fixture> --kernel-blob <sha>\n", argv[0]);
    return 2;
}
int main(int argc, char** argv) {
    try { return run(argc, argv); }
    catch (const Refusal& e) { std::fflush(stdout); std::fprintf(stderr, "REPLAY ERROR: %s\n", e.msg.c_str()); return 2; }
}
