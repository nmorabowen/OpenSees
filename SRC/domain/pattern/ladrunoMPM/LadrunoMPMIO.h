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

// Authors: Nicolas Mora Bowen, Guppi (Ladruno)
// Created: 11/2026
//
// LadrunoMPMIO — the two MPM file contracts on the fork side (WP-232c phase c2;
// Ladruno_implementation/232a_mpm_contracts.md):
//
//   * readParticleSet  — Contract 1, the particle set apeGmsh writes
//     (`<deck stem>.particles.h5`), read by `mpm particles -file <path>`. Every
//     rule of the contract is checked before anything is loaded, and every
//     problem is reported with its rule code (C1-...), as the executable spec
//     tests/mpm_contracts.py::validate_particle_set does.
//   * ParticleWriter   — Contract 2, the particle group of a `.ladruno` file
//     (FORMAT_VERSION 1, FORMAT_MINOR_VERSION 2), driven by `mpm output`. Its
//     stages are particle-only, written exactly as LadrunoRecorder writes an
//     empty stage plus MODEL/PARTICLES and RESULTS/ON_PARTICLES. Checked by
//     Ladruno_scripts/ladruno_recorder_tests/ladruno_format.py::validate.
//
// Both need HDF5 (`_HDF5`, the LadrunoRecorder build condition); without it
// they refuse with a message.

#ifndef LadrunoMPMIO_h
#define LadrunoMPMIO_h

#include <string>
#include <vector>

namespace ladruno_mpm_io {

bool hdf5Available(void);

// ---- Contract 1 -------------------------------------------------------------
struct ParticleSet {
  std::string path;                     // as opened
  std::string setId, generator;
  std::string unitLength, unitForce, unitTime, unitSystem;   // unitSystem may be ""
  int    ndm = 0;
  long long n = 0;
  double thickness = 0.0;               // 2D only
  double origin[3] = {0, 0, 0}, cellSize[3] = {0, 0, 0};
  int    numCells[3] = {0, 0, 0};
  std::string basis;
  std::vector<int> bodyId;
  std::vector<std::string> bodyName;
  std::vector<long long> id;
  std::vector<int> body, mat;
  std::vector<double> x, volume, mass, lp, v, stress, stress33;   // mass = MASS or DENSITY*VOLUME
  bool hasVelocity = false, hasStress = false, hasStress33 = false;
};

// 0 when the file honours every rule; otherwise -1 and `problems` lists every
// broken rule, each message starting with its code.
int readParticleSet(const std::string &path, ParticleSet &out, std::vector<std::string> &problems);

// ---- Contract 2 -------------------------------------------------------------
struct WriterConfig {
  std::string path;
  int    every = 1;          // output every N steps (when dtOut <= 0)
  double dtOut = -1.0;       // or every dtOut of Domain time
  int    compress = 1;       // deflate level (0: no filter)
  double flushSeconds = 0.0;    // wall-clock H5Fflush cadence (LadrunoRecorder -flush, WP-164);
                                // <= 0 (the default): after every output -- particle outputs are
                                // already throttled by -every / -dt, and a file is crash-
                                // consistent only at a flush (c2 review 2)
  std::vector<std::string> fields;   // POSITION always first
};

struct StageMeta {
  int    epoch = 0;
  int    patternTag = -1;
  std::string mode = "explicit";
  double thickness = 1.0;   // written when ndm = 2
  std::string unitLength, unitForce, unitTime, unitSystem;
  std::vector<std::string> sourceSetIds;
  double origin[3] = {0, 0, 0}, cellSize[3] = {0, 0, 0};
  int    numCells[3] = {0, 0, 0};
  std::string basis;
  std::vector<int> bodyId;
  std::vector<std::string> bodyName;
  // the registry, sorted by id
  std::vector<long long> id;
  std::vector<int> body, mat;
  std::vector<double> mass;
};

class ParticleWriter {
 public:
  ParticleWriter();
  ~ParticleWriter();
  // Configure the writer. The file (and its INFO group) is created only when
  // the first stage opens, so a writer that never outputs leaves no file
  // (an INFO-only file would not validate). -1 with `err` on bad input.
  int  open(const WriterConfig &cfg, int ndm, std::string &err);
  bool isConfigured(void) const;
  bool isOpen(void) const;      // the file exists
  // Open MODEL_STAGE[<next>] at commit tag `step`, time `time`.
  int  openStage(const StageMeta &meta, int step, double time, std::string &err);
  bool stageOpen(void) const;
  int  stageEpoch(void) const;
  // One output: ids (strictly increasing), and for each configured field its
  // rows (ids.size() x nComp, row-major, in the order of cfg.fields).
  int  writeOutput(int step, double time, const std::vector<long long> &ids,
                   const std::vector<std::vector<double> > &fieldRows, std::string &err);
  int  numOutputs(void) const;
  void flush(void);
  // Flush when cfg.flushSeconds of wall clock passed since the last flush: a
  // killed run then leaves a readable file (an HDF5 file never flushed is not).
  void flushIfDue(void);
  void close(void);
  const WriterConfig &config(void) const;
  // Components of a reserved field for ndm (232a §2.6), or "" if not written by v1.
  static std::string components(const std::string &field, int ndm);

  struct Impl;              // the HDF5 state (LadrunoMPMIO.cpp)

 private:
  Impl *impl;
};

}  // namespace ladruno_mpm_io

#endif
