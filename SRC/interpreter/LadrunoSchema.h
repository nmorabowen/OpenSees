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

// Ladruno WP-185 (apeGmsh#1667, program link F3): `ladrunoSchema` -- a
// machine-readable, deterministic description of what this build offers, so a
// frontend (apeGmsh's BackendInfo) can detect capabilities instead of keying on
// the build sha alone.
//
// One JSON object, keys sorted at every level, compact separators, no paths and
// no timestamps (byte-identical across runs of one binary):
//
//   {"build":"<ladrunoBuild>",
//    "capabilities":["mdStageFlagPreserved",...],   sorted; kCapabilities below
//    "commands":["LadrunoStaggeredAnalyze",...],   sorted; the rows of
//                                                   LadrunoCommandTable.h this
//                                                   ENGINE registers
//    "engine":"dl"|"classic",                       dl = TclWrapper + Python,
//                                                   classic = OpenSees.exe/SP/MP
//    "options":{"LadrunoRCConcrete":{"-Cd":{"nargs":"list","since":"PR-155"},..},..},
//    "schema_version":1}
//
// Everything is GENERATED from the declarations: `commands` from the command
// table (X-macro, stringified columns, so no entry point is referenced and the
// header is safe in both the classic and the DL translation units), `options`
// from ladruno_opt::kOptTables (SRC/material/LadrunoOptSpec.h), `capabilities`
// from kCapabilities below. `nargs` is the declared value count, or "list" for a
// variable-length list.
//
// Consumers: a token is supported iff it is present in options[<owner>]; a
// command is registered by this engine iff it is in `commands`; a behaviour
// holds iff its name is in `capabilities`. `since` is informational provenance
// (PR-155, ADR19-2b.2c.1, ...), not an orderable version; never compare two
// `since` values.
//
// Versions (apeGmsh#2517): v0 had exactly the keys build, commands, engine,
// options, schema_version, and both its readers (tests/test_wp185_ladruno_schema.py,
// apeGmsh BackendInfo) refuse any other key set. v1 adds `capabilities`, so it is
// a new version, not a silent v0 extension. From v1 on: a new ROW in an existing
// section (a command, an option token, a capability) does not bump the version;
// a new top-level section or a changed value shape does.
//
// Fail closed: a duplicate command, capability or (owner, token) key would make
// the JSON ambiguous, so render() refuses instead of emitting it. Arguments to
// the command itself are refused by the two bridges (OpenSeesOutputCommands.cpp,
// SRC/tcl/commands.cpp), which are the only includers.

#ifndef LadrunoSchema_h
#define LadrunoSchema_h

#include <LadrunoOptSpec.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <string>
#include <vector>

namespace ladruno_schema {

enum Engine { ENGINE_DL, ENGINE_CLASSIC };

static const int kSchemaVersion = 1;

// Behavioural capabilities: facts about this build that no command and no option
// token shows, so a consumer cannot otherwise tell a build that has a fix from one
// that does not (a bare `ladrunoBuild` sha cannot be ordered against a commit).
// A name here is a PROMISE about the code: add it in the commit that makes it true
// (or after), never before, and give it a test that fails when the behaviour
// breaks while the name is still listed. Never remove a name whose behaviour still
// holds; consumers read absence as "older build, refuse".
static const char* const kCapabilities[] = {
  // ManzariDafalias' class-wide stage flag (static mElastFlag) survives the
  // construction of any further ManzariDafalias or LadrunoSANISAND instance:
  // no constructor (full, or the classTag one getCopy() runs) resets it, so an
  // element activation or a new nDMaterial after `updateMaterialStage ... -stage 1`
  // keeps every MD material plastic. Fix 43b235bbb (fork PR #1005, apeGmsh#2444);
  // marker apeGmsh#2517. Test: tests/test_md_stage_flag_capability.py.
  "mdStageFlagPreserved",
};

inline void appendJsonString(std::string& o, const char* s)
{
  o += '"';
  for (const char* p = s ? s : ""; *p; ++p) {
    const unsigned char c = (unsigned char)*p;
    if (c == '"' || c == '\\') { o += '\\'; o += (char)c; }
    else if (c < 0x20) {
      char buf[8];
      snprintf(buf, sizeof(buf), "\\u%04x", (unsigned)c);
      o += buf;
    }
    else o += (char)c;
  }
  o += '"';
}

inline bool lessStr(const char* a, const char* b) { return strcmp(a, b) < 0; }

// Build the schema for `engine`. Returns false (and sets `err`) on a duplicate
// key, never a half-written or ambiguous document.
inline bool render(Engine engine, std::string& out, std::string& err)
{
  struct Row { const char* name; const char* dl; const char* classic; };
  static const Row rows[] = {
#define LADRUNO_COMMAND(name, dl, classic) {name, #dl, #classic},
#include "LadrunoCommandTable.h"
#undef LADRUNO_COMMAND
  };

  std::vector<const char*> cmds;
  for (const Row& r : rows) {
    const char* col = (engine == ENGINE_CLASSIC) ? r.classic : r.dl;
    if (strcmp(col, "LADRUNO_NONE") != 0) cmds.push_back(r.name);
  }
  std::sort(cmds.begin(), cmds.end(), lessStr);
  for (size_t i = 1; i < cmds.size(); ++i)
    if (strcmp(cmds[i - 1], cmds[i]) == 0) {
      err = std::string("duplicate command '") + cmds[i] + "' in LadrunoCommandTable.h";
      return false;
    }

  std::vector<const char*> caps(std::begin(kCapabilities), std::end(kCapabilities));
  std::sort(caps.begin(), caps.end(), lessStr);
  for (size_t i = 1; i < caps.size(); ++i)
    if (strcmp(caps[i - 1], caps[i]) == 0) {
      err = std::string("duplicate capability '") + caps[i] + "' in ladruno_schema::kCapabilities";
      return false;
    }

  const size_t nTables = sizeof(ladruno_opt::kOptTables) / sizeof(ladruno_opt::kOptTables[0]);
  std::vector<const ladruno_opt::OptTable*> tables;
  for (size_t i = 0; i < nTables; ++i) tables.push_back(&ladruno_opt::kOptTables[i]);
  std::sort(tables.begin(), tables.end(),
            [](const ladruno_opt::OptTable* a, const ladruno_opt::OptTable* b) {
              return strcmp(a->owner, b->owner) < 0; });
  for (size_t i = 1; i < tables.size(); ++i)
    if (strcmp(tables[i - 1]->owner, tables[i]->owner) == 0) {
      err = std::string("duplicate owner '") + tables[i]->owner + "' in ladruno_opt::kOptTables";
      return false;
    }

  std::string o;
  o += "{\"build\":";
#ifdef OPENSEES_VERSION
  appendJsonString(o, OPENSEES_VERSION);
#else
  appendJsonString(o, "unknown");
#endif
  o += ",\"capabilities\":[";
  for (size_t i = 0; i < caps.size(); ++i) {
    if (i) o += ',';
    appendJsonString(o, caps[i]);
  }
  o += "],\"commands\":[";
  for (size_t i = 0; i < cmds.size(); ++i) {
    if (i) o += ',';
    appendJsonString(o, cmds[i]);
  }
  o += "],\"engine\":";
  appendJsonString(o, engine == ENGINE_CLASSIC ? "classic" : "dl");
  o += ",\"options\":{";
  for (size_t t = 0; t < tables.size(); ++t) {
    const ladruno_opt::OptTable& tab = *tables[t];
    std::vector<const ladruno_opt::OptSpec*> specs;
    for (size_t k = 0; k < tab.n; ++k) specs.push_back(&tab.opts[k]);
    std::sort(specs.begin(), specs.end(),
              [](const ladruno_opt::OptSpec* a, const ladruno_opt::OptSpec* b) {
                return strcmp(a->name, b->name) < 0; });
    for (size_t k = 1; k < specs.size(); ++k)
      if (strcmp(specs[k - 1]->name, specs[k]->name) == 0) {
        err = std::string("duplicate option '") + specs[k]->name + "' declared for " + tab.owner;
        return false;
      }
    if (t) o += ',';
    appendJsonString(o, tab.owner);
    o += ":{";
    for (size_t k = 0; k < specs.size(); ++k) {
      if (k) o += ',';
      appendJsonString(o, specs[k]->name);
      o += ":{\"nargs\":";
      if (specs[k]->nargs == ladruno_opt::LIST) o += "\"list\"";
      else o += std::to_string(specs[k]->nargs);
      o += ",\"since\":";
      appendJsonString(o, specs[k]->since);
      o += '}';
    }
    o += '}';
  }
  o += "},\"schema_version\":";
  o += std::to_string(kSchemaVersion);
  o += '}';

  out.swap(o);
  return true;
}

}  // namespace ladruno_schema

#endif
