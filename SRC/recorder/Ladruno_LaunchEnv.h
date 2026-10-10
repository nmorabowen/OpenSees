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

#ifndef Ladruno_LaunchEnv_h
#define Ladruno_LaunchEnv_h

/*************************************************************************************

Ladruno_LaunchEnv.h  (WP-163 M4 / MP-3)

Which rank of an MPI launch is this process? The recorders live in the shared
sequential OPS_Recorder library (no MPI — LEDGER_quirks "OPS_Recorder is compiled
sequentially"), so on the interpreter-per-rank path (openseesmp) they read the
launcher's per-rank environment. One probe, shared by LadrunoRecorder and the
EnergyBalance recorder (they had diverging copies).

Pairs, first with SIZE > 1 wins:
  PMI_SIZE / PMI_RANK                          Intel MPI, MS-MPI, MPICH hydra
  OMPI_COMM_WORLD_SIZE / OMPI_COMM_WORLD_RANK  OpenMPI
  SLURM_NTASKS / SLURM_PROCID                  srun (any MPI plugin, incl. pmix)

Two guards the old copies lacked:
  * SLURM is trusted ONLY inside an srun job step. `sbatch --ntasks=N` exports
    SLURM_NTASKS=N (and SLURM_PROCID=0) into the batch shell itself, so a plain
    sequential run in a batch script was written as `part-0` of N. srun sets
    SLURM_STEP_ID to a real step number; the batch/extern pseudo-steps use the
    reserved values 0xFFFFFFFE / 0xFFFFFFFD (or leave it unset).
  * A SIZE > 1 with a missing or out-of-range RANK is an ERROR, not rank 0:
    every rank would otherwise truncate the same `part-0` file.

**************************************************************************************/

#include <cstdlib>
#include <cerrno>
#include <string>

namespace ladruno {
namespace launch {

	enum Status { NotLaunched = 0, Launched = 1, Inconsistent = -1 };

	// Parse a non-negative decimal integer; false on anything else.
	inline bool parseNonNegative(const char* s, long long& out)
	{
		if (s == 0 || *s == '\0')
			return false;
		char* end = 0;
		errno = 0;
		long long v = std::strtoll(s, &end, 10);
		if (errno != 0 || end == s || *end != '\0' || v < 0)
			return false;
		out = v;
		return true;
	}

	// NotLaunched: single process (no launcher env, or SIZE <= 1).
	// Launched:    rank in [0, size), size > 1; `source` names the env pair.
	// Inconsistent: a launcher claims SIZE > 1 but RANK is missing/invalid;
	//               `error` says which variable — callers must not guess rank 0.
	inline Status detectRank(int& rank, int& size, std::string& source, std::string& error)
	{
		static const char* const pairs[][2] = {
			{ "PMI_SIZE",             "PMI_RANK" },
			{ "OMPI_COMM_WORLD_SIZE", "OMPI_COMM_WORLD_RANK" },
			{ "SLURM_NTASKS",         "SLURM_PROCID" },
		};
		rank = 0;
		size = 1;
		for (size_t i = 0; i < sizeof(pairs) / sizeof(pairs[0]); ++i) {
			long long np = 0;
			if (!parseNonNegative(std::getenv(pairs[i][0]), np) || np <= 1)
				continue;
			if (i == 2) {
				// SLURM: only inside a real srun step (see the header note).
				long long step = 0;
				if (!parseNonNegative(std::getenv("SLURM_STEP_ID"), step) ||
				    step >= 4294967290LL)
					continue;
			}
			long long r = 0;
			if (!parseNonNegative(std::getenv(pairs[i][1]), r) || r >= np) {
				error = std::string(pairs[i][0]) + " > 1 but " + pairs[i][1] +
				        " is missing or not in [0, " + pairs[i][0] + ")";
				return Inconsistent;
			}
			rank = (int)r;
			size = (int)np;
			source = std::string(pairs[i][0]) + "/" + pairs[i][1];
			return Launched;
		}
		return NotLaunched;
	}

} // namespace launch
} // namespace ladruno

#endif // Ladruno_LaunchEnv_h
