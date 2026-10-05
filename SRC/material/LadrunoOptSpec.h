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

// Ladruno (WP-167): the fail-closed unknown-token policy for fork parsers.
//
// WHY. A fork parser whose option ladder simply falls through on a token it
// does not recognise (the old LadrunoRCConcrete loop ended in a comment calling
// that "forward-compat") turns a misspelled flag, or a flag this build does not
// have yet, into a SILENT no-op: the model runs with the default physics and
// nothing says so. apeGmsh #1184 emitted `-crackedNu`/`-betaC` for
// LadrunoRCConcrete; this build dropped both without a word. The policy here is
// the opposite:
//
//   * every option a parser accepts is DECLARED in a table: its name, how many
//     values follow it, and `since` -- the fork PR/WP/phase that introduced it,
//     so a refusal can tell the user which build they need;
//   * a token that is not in the table makes the command FAIL (the parser
//     returns null, so Tcl raises an error and openseespy an exception) and the
//     message names the token, the material, its tag and this build's stamp.
//
// A frontend that wants to emit a NEW option must wait for the build that
// declares it; feature detection is `ladrunoBuild` + the `since` field, never
// "send it and hope it is ignored".
//
// Usage (one table per parser family; the ladder that consumes the values
// stays as it was, the table only GATES it):
//
//   static const ladruno_opt::OptSpec kOpts[] = { {"-rho", 1, "PR-155"}, ... };
//   const char* opt = OPS_GetString();
//   const ladruno_opt::OptSpec* s = ladruno_opt::find(kOpts, opt);
//   if (s == 0)  { ladruno_opt::reportUnknown("MyMat", tag, opt, kOpts); return 0; }
//   if (!ladruno_opt::haveValues("MyMat", tag, *s)) return 0;
//   if (strcmp(opt, "-rho") == 0) ...      // the existing ladder
//   else { ladruno_opt::reportUnhandled("MyMat", tag, opt); return 0; }
//
// The final `else` makes the table and the ladder fail closed in BOTH
// directions: a ladder branch missing from the table is refused as unknown
// (the parser's enumerating test catches it), and a table entry with no branch
// is refused as unhandled instead of being skipped.
//
// The quirk lint (`ci/check_quirk_patterns.py`, rule `unknown-token`) flags, in
// fork parsers, a comment that says unrecognised tokens are let through and an
// option ladder that is the last statement of its loop with no final `else`.

#ifndef LadrunoOptSpec_h
#define LadrunoOptSpec_h

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <stddef.h>
#include <string.h>

namespace ladruno_opt {

// nargs: number of values that follow the option (0 = a bare flag), or
// LIST for a variable-length numeric list that the parser reads until the
// next non-number (e.g. a backbone `-Ce e1 e2 ...`).
enum { LIST = -1 };

struct OptSpec {
  const char* name;    // the token, with its leading '-'
  int         nargs;   // values that follow it, or LIST
  const char* since;   // the fork PR/WP/phase that introduced it
};

inline const char* buildStamp()
{
#ifdef OPENSEES_VERSION
  return OPENSEES_VERSION;
#else
  return "unknown";
#endif
}

template <size_t N>
inline const OptSpec* find(const OptSpec (&table)[N], const char* tok)
{
  if (tok == 0) return 0;
  for (size_t i = 0; i < N; ++i)
    if (strcmp(table[i].name, tok) == 0) return &table[i];
  return 0;
}

// The refusal. The first line names the token (what a caller greps for) and the
// build; the next lists the accepted options with their `since`.
template <size_t N>
inline void reportUnknown(const char* material, int tag, const char* tok,
                          const OptSpec (&table)[N])
{
  opserr << "ERROR nDMaterial " << material << " " << tag << ": unknown option '"
         << (tok ? tok : "(null)") << "' (build " << buildStamp() << ")\n"
         << "  unknown options are fatal (WP-167): an option this build does not declare is "
            "refused, never ignored.\n"
         << "  accepted:";
  for (size_t i = 0; i < N; ++i) {
    opserr << " " << table[i].name;
    if (table[i].nargs == LIST)  opserr << " {..}";
    else if (table[i].nargs > 0) opserr << " (" << table[i].nargs << ")";
    opserr << " [" << table[i].since << "]";
  }
  opserr << "\n";
}

// The option is declared but too few tokens remain for its values.
inline bool haveValues(const char* material, int tag, const OptSpec& s)
{
  if (s.nargs <= 0 || OPS_GetNumRemainingInputArgs() >= s.nargs) return true;
  opserr << "ERROR nDMaterial " << material << " " << tag << ": option '" << s.name
         << "' needs " << s.nargs << " value(s), " << OPS_GetNumRemainingInputArgs()
         << " left (build " << buildStamp() << ")\n";
  return false;
}

// A declared option the parser's ladder has no branch for (a programming error:
// the table and the ladder drifted). Refused rather than skipped.
inline void reportUnhandled(const char* material, int tag, const char* tok)
{
  opserr << "ERROR nDMaterial " << material << " " << tag << ": option '" << tok
         << "' is declared but has no handler in this build (build " << buildStamp()
         << ") -- the option table and the parser disagree\n";
}

// ---------------------------------------------------------------------------
// The LadrunoRCConcrete family (LadrunoRCConcrete 33015, LadrunoRCFiniteStrain
// 33018): ONE grammar, so ONE table. `since` = the PR / ADR-19 phase that added
// the option to the family. -crackedNu and -betaC (apeGmsh #1184, open fork PR
// #877) are deliberately ABSENT: they are refused until that PR declares them
// here with its own `since`.
// ---------------------------------------------------------------------------
static const OptSpec kLadrunoRCOptions[] = {
  // Phase 1 (#155/#192): backbones, compression softening, tangent mode
  {"-Ce", LIST, "PR-155"}, {"-Cs", LIST, "PR-155"}, {"-Cd", LIST, "PR-155"},
  {"-Te", LIST, "PR-155"}, {"-Ts", LIST, "PR-155"}, {"-Td", LIST, "PR-155"},
  {"-Kc", 1, "PR-155"}, {"-betaFloor", 1, "PR-155"}, {"-rho", 1, "PR-155"},
  {"-beta", 0, "PR-155"}, {"-lublinerReduced", 0, "PR-155"},
  {"-secant", 0, "PR-192"}, {"-numericalTangent", 0, "PR-192"},
  // Phase 2a (#239): fixed-crack aggregate interlock
  {"-interlock", 0, "PR-239"}, {"-agg", 1, "PR-239"}, {"-crackStrain", 1, "PR-239"},
  {"-crackSpacing", 1, "PR-239"}, {"-lch", 1, "PR-239"}, {"-betaSrMin", 1, "PR-239"},
  // Phase 2b.1 (#245): cyclic friction-slip
  {"-cyclic", 0, "PR-245"},
  // Phase 2b.2b (#253): X-cracking + interlock wear
  {"-xcrack", 0, "PR-253"}, {"-degKappa", 1, "PR-253"}, {"-degSlipRef", 1, "PR-253"},
  {"-degMin", 1, "PR-253"},
  // Phase 4a (#263): IMPL-EX
  {"-implex", 0, "PR-263"}, {"-implexAlpha", 1, "PR-263"}, {"-implexControl", 2, "PR-263"},
  // Phase 2b.2c.1 (62067bae3): crack-shear retention curves
  {"-shearRetention", 1, "ADR19-2b.2c.1"}, {"-shearRetFactor", 1, "ADR19-2b.2c.1"},
  // Phase 3a (#273): tension stiffening
  {"-tensStiff", 1, "PR-273"}, {"-tensStiffC", 1, "PR-273"}, {"-tensStiffAlpha", 1, "PR-273"},
  // Phase 3b (#277): crack-band regularization
  {"-autoRegularization", 1, "PR-277"},
};

}  // namespace ladruno_opt

#endif
