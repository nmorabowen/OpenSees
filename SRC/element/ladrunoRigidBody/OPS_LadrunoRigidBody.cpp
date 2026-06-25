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

// Ladruno: OPS parser for the RigidBody element (ADR 58).
//
//   element LadrunoRigidBody tag  N  s1 ... sN   [-mass m] [-internalNode tag]
//
//   N           number of slave nodes that form the rigid body
//   s1..sN      the slave node tags (3- or 6-DOF, 3D)
//   -mass m     total body mass (else condensed from the slaves' nodal mass)
//   -internalNode tag   fix the private CoM node tag (else auto-assigned)
//
//   The body moves as one 6-DOF rigid object; its mass is condensed to a private
//   internal node at the center of mass and the slaves become rigid followers.
//   P1 = translation/ballistic (R=I); SO(3) rotation is added at P2. See ADR 58.
//
// Written: N. Mora-Bowen (Ladruno), 2026.

#include <LadrunoRigidBody.h>
#include <ID.h>
#include <Domain.h>
#include <Element.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>

void* OPS_LadrunoRigidBody(void)
{
  int ndm = OPS_GetNDM();
  if (ndm != 3) {
    opserr << "WARNING LadrunoRigidBody: model ndm must be 3 (v1, P1)\n";
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING insufficient args\n"
           << "Want: element LadrunoRigidBody tag N s1..sN [-mass m] "
           << "[-internalNode tag]\n";
    return 0;
  }

  int idata[2];                    // tag, N
  int n = 2;
  if (OPS_GetIntInput(&n, idata) < 0) {
    opserr << "WARNING LadrunoRigidBody: invalid tag/N\n";
    return 0;
  }
  int tag = idata[0], N = idata[1];
  if (N < 1) {
    opserr << "WARNING LadrunoRigidBody: N (number of slave nodes) must be >= 1; got "
           << N << "\n";
    return 0;
  }
  if (OPS_GetNumRemainingInputArgs() < N) {
    opserr << "WARNING LadrunoRigidBody: need " << N << " slave node tags\n";
    return 0;
  }
  ID slaves(N);
  for (int i = 0; i < N; i++) {
    int h; n = 1;
    if (OPS_GetIntInput(&n, &h) < 0) {
      opserr << "WARNING LadrunoRigidBody: invalid slave node " << i << "\n";
      return 0;
    }
    slaves(i) = h;
  }

  double mUser = -1.0;
  int intNodeTag = -1;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char* opt = OPS_GetString();
    if (strcmp(opt, "-mass") == 0) {
      n = 1;
      if (OPS_GetDoubleInput(&n, &mUser) < 0) {
        opserr << "WARNING LadrunoRigidBody: -mass wants a value\n";
        return 0;
      }
      if (mUser < 0.0) {
        opserr << "WARNING LadrunoRigidBody: -mass must be >= 0\n";
        return 0;
      }
    }
    else if (strcmp(opt, "-internalNode") == 0 || strcmp(opt, "-intNode") == 0) {
      n = 1;
      if (OPS_GetIntInput(&n, &intNodeTag) < 0) {
        opserr << "WARNING LadrunoRigidBody: -internalNode wants a node tag\n";
        return 0;
      }
    }
    else {
      opserr << "WARNING LadrunoRigidBody: unknown option '" << opt << "'\n";
      return 0;
    }
  }

  Element* e = new LadrunoRigidBody(tag, ndm, slaves, mUser, intNodeTag);
  if (e == 0) {
    opserr << "WARNING LadrunoRigidBody: could not create element\n";
    return 0;
  }
  return e;
}
