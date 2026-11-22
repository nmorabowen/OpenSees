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

// Implementation of LadrunoNorSandPlaneStrain; see LadrunoNorSandPlaneStrain.h and LadrunoNorSand.h.
// Written: N. Mora-Bowen (Ladruno), 2026.

#include "LadrunoNorSandPlaneStrain.h"

LadrunoNorSandPlaneStrain::LadrunoNorSandPlaneStrain()
  : LadrunoNorSand(ND_TAG_LadrunoNorSandPlaneStrain, LadrunoNorSand::DIM_PSTRAIN)
{
}

LadrunoNorSandPlaneStrain::LadrunoNorSandPlaneStrain(int tag, const ladruno_norsand::Params& p, const double sig0[6],
                                       double v0, double pi0, double dens)
  : LadrunoNorSand(tag, ND_TAG_LadrunoNorSandPlaneStrain, p, sig0, v0, pi0, dens, LadrunoNorSand::DIM_PSTRAIN)
{
}

LadrunoNorSandPlaneStrain::~LadrunoNorSandPlaneStrain()
{
}

NDMaterial*
LadrunoNorSandPlaneStrain::getCopy(void)
{
  LadrunoNorSandPlaneStrain* clone = new LadrunoNorSandPlaneStrain();
  clone->copyFrom(*this);
  return clone;
}
