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

// Ladruno WP-185 (apeGmsh#1667, program link F3): `ladrunoSchema` v0 -- a
// machine-readable, deterministic description of what this build offers, so a
// frontend (apeGmsh's BackendInfo) can detect capabilities instead of keying on
// the build sha alone.
//
// One JSON object, keys sorted at every level, compact separators, no paths and
// no timestamps (byte-identical across runs of one binary):
//
//   {"build":"<ladrunoBuild>",
//    "commands":["LadrunoStaggeredAnalyze",...],   sorted; the rows of
//                                                   LadrunoCommandTable.h this
//                                                   ENGINE registers
//    "engine":"dl"|"classic",                       dl = TclWrapper + Python,
//                                                   classic = OpenSees.exe/SP/MP
//    "options":{"LadrunoRCConcrete":{"-Cd":{"nargs":"list","since":"PR-155"},..},..},
//    "schema_version":0}
//
// Everything is GENERATED from the declarations: `commands` from the command
// table (X-macro, stringified columns, so no entry point is referenced and the
// header is safe in both the classic and the DL translation units), `options`
// from ladruno_opt::kOptTables (SRC/material/LadrunoOptSpec.h). `nargs` is the
// declared value count, or "list" for a variable-length list.
//
// Consumers: a token is supported iff it is present in options[<owner>].
// `since` is informational provenance (PR-155, ADR19-2b.2c.1, ...), not an
// orderable version; never compare two `since` values.
//
// Fail closed: a duplicate command or a duplicate (owner, token) key would make
// the JSON ambiguous, so render() refuses instead of emitting it. Arguments to
// the command itself are refused by the two bridges (OpenSeesOutputCommands.cpp,
// SRC/tcl/commands.cpp), which are the only includers.

#ifndef LadrunoSchema_h
#define LadrunoSchema_h

#include <LadrunoOptSpec.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace ladruno_schema {

enum Engine { ENGINE_DL, ENGINE_CLASSIC };

static const int kSchemaVersion = 0;

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
  o += ",\"commands\":[";
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
