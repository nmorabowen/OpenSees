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

// Unified 2D/3D dispatcher for `element LadrunoForceBeamColumn`.
// Inspects OPS_GetNDM() / OPS_GetNDF() and forwards to the 2D or 3D
// parser in LadrunoForceBeamColumn2d.cpp / LadrunoForceBeamColumn3d.cpp.

#include <elementAPI.h>
#include <OPS_Globals.h>

extern void* OPS_LadrunoForceBeamColumn2d();
extern void* OPS_LadrunoForceBeamColumn3d();

void* OPS_LadrunoForceBeamColumn()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if (ndm == 2 && ndf == 3) {
        return OPS_LadrunoForceBeamColumn2d();
    } else if (ndm == 3 && ndf == 6) {
        return OPS_LadrunoForceBeamColumn3d();
    } else {
        opserr << "WARNING LadrunoForceBeamColumn -- model must be ndm 2/ndf 3 (2D) "
               << "or ndm 3/ndf 6 (3D); got ndm=" << ndm << " ndf=" << ndf << "\n";
        return 0;
    }
}
