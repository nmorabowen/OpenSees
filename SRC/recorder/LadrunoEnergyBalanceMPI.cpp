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

// Ladruno -- see LadrunoEnergyBalanceMPI.h for why this is its own TU.

#include "LadrunoEnergyBalanceMPI.h"

#if defined(_PARALLEL_INTERPRETERS)
#include <mpi.h>

int
ladrunoEnergyBalanceAgreeMax(int *flags, int n)
{
  if (n <= 0) return 0;
  int init = 0, fin = 0;
  MPI_Initialized(&init);
  if (!init) return 0;
  MPI_Finalized(&fin);
  if (fin) return 0;
  int np = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  if (np <= 1) return 0;
  int local[16];
  int m = (n > 16) ? 16 : n;
  for (int i = 0; i < m; i++) local[i] = flags[i];
  MPI_Allreduce(local, flags, m, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return 0;
}

#else

int
ladrunoEnergyBalanceAgreeMax(int *, int)
{
  return 0;
}

#endif
