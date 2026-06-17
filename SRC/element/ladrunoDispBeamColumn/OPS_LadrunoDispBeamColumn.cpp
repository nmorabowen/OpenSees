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

// Author: N. Mora-Bowen (Ladruno), 06/2026
//
// ndm-dispatch factory for the LadrunoDispBeamColumn element (Tcl + Python).
// One user-facing command picks the planar or spatial class from the model
// dimension (same UX as dispBeamColumn / elasticBeamColumn). See
// Ladruno_implementation/32_ladruno_dispbeamcolumn_regularization_adr.md.
//
// Usage (2D):
//   element('LadrunoDispBeamColumn', eleTag, iNode, jNode, transfTag, integrationTag
//           [, '-mass', massPerLength] [, '-cMass'] [, '-damp', dampTag]
//           [, '-lch', {ip|element|<value>}])   # regularization length (default: ip)
//
// -lch selects the characteristic length reported to crack-band / auto-regularizing
// materials (ASDConcrete1D -autoRegularization, ASDSteel1D, LadrunoUniaxialJ2+Lemaitre):
//   ip       per-integration-point tributary length wt[i]*L (default, mesh-objective)
//   element  the whole element length L (A/B-debug only; warns about energy double-count)
//   <value>  a fixed user crack-band width
//
// NOTE (Stage 1): the 3D sibling (LadrunoDispBeamColumn3d) is not yet built; a 3D
// model currently errors here. Tier-2 (embedded softening hinge) is a later stage.

#include <elementAPI.h>
#include <OPS_Globals.h>

// Planar builder, defined in LadrunoDispBeamColumn2d.cpp.
extern void *OPS_LadrunoDispBeamColumn2d();

void *OPS_LadrunoDispBeamColumn()
{
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();

  if (ndm == 2 && ndf == 3)
    return OPS_LadrunoDispBeamColumn2d();

  opserr << "WARNING LadrunoDispBeamColumn -- model must be ndm 2 / ndf 3 (2D). "
            "The 3D sibling is not yet available in this build.\n";
  return 0;
}
