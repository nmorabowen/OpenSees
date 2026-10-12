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

// Ladruno WP-188 -- see LadrunoFactorReuseMPI.h for why this is its own TU.

#include "LadrunoFactorReuseMPI.h"

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
#include <mpi.h>

static bool ladrunoMpiLive(void)
{
  int init = 0, fin = 0;
  MPI_Initialized(&init);
  if (!init) return false;
  MPI_Finalized(&fin);
  return fin == 0;
}

int
ladrunoFactorReuseAgreeMax(int *flags, int n)
{
  if (n <= 0 || !ladrunoMpiLive()) return 0;
  int np = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  if (np <= 1) return 0;
  int local[8];
  int m = (n > 8) ? 8 : n;
  for (int i = 0; i < m; i++) local[i] = flags[i];
  MPI_Allreduce(local, flags, m, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return 0;
}

int
ladrunoFactorReuseRank(void)
{
  if (!ladrunoMpiLive()) return 0;
  int r = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &r);
  return r;
}

#else

int
ladrunoFactorReuseAgreeMax(int *, int)
{
  return 0;
}

int
ladrunoFactorReuseRank(void)
{
  return 0;
}

#endif
