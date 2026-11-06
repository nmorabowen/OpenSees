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

// Implementation of LadrunoNorSand3D; see LadrunoNorSand3D.h and LadrunoNorSand.h.
// Written: N. Mora-Bowen (Ladruno), 2026.

#include "LadrunoNorSand3D.h"

LadrunoNorSand3D::LadrunoNorSand3D()
  : LadrunoNorSand(ND_TAG_LadrunoNorSand3D, LadrunoNorSand::DIM_3D)
{
}

LadrunoNorSand3D::LadrunoNorSand3D(int tag, const ladruno_norsand::Params& p, const double sig0[6],
                                       double v0, double pi0, double dens)
  : LadrunoNorSand(tag, ND_TAG_LadrunoNorSand3D, p, sig0, v0, pi0, dens, LadrunoNorSand::DIM_3D)
{
}

LadrunoNorSand3D::~LadrunoNorSand3D()
{
}

NDMaterial*
LadrunoNorSand3D::getCopy(void)
{
  LadrunoNorSand3D* clone = new LadrunoNorSand3D();
  clone->copyFrom(*this);
  return clone;
}
