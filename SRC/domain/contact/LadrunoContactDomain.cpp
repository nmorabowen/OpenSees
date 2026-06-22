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

// ADR-39 P1b — LadrunoContactDomain implementation. See header for design.

#include "LadrunoContactDomain.h"
#include "LadrunoContactSurface.h"
#include <OPS_Globals.h>   // opserr, endln

LadrunoContactDomain::LadrunoContactDomain()
  : numCommits(0), numReverts(0)
{
}

LadrunoContactDomain::~LadrunoContactDomain()
{
    for (size_t i = 0; i < theSurfaces.size(); i++)
        delete theSurfaces[i];
    theSurfaces.clear();
    theContacts.clear();
}

int
LadrunoContactDomain::addSurface(LadrunoContactSurface *surf)
{
    if (surf == 0)
        return -1;
    // reject a duplicate tag
    if (getSurface(surf->getTag()) != 0) {
        opserr << "WARNING LadrunoContactDomain::addSurface() - surface tag "
               << surf->getTag() << " already exists\n";
        return -1;
    }
    theSurfaces.push_back(surf);
    return 0;
}

LadrunoContactSurface *
LadrunoContactDomain::getSurface(int tag) const
{
    for (size_t i = 0; i < theSurfaces.size(); i++)
        if (theSurfaces[i]->getTag() == tag)
            return theSurfaces[i];
    return 0;
}

int
LadrunoContactDomain::addContact(int tag, int masterSurfTag, int slaveSurfTag,
                                 double kn, double kt, double mu)
{
    if (getSurface(masterSurfTag) == 0 || getSurface(slaveSurfTag) == 0) {
        opserr << "WARNING LadrunoContactDomain::addContact() - master/slave surface "
                  "not defined (master " << masterSurfTag << ", slave " << slaveSurfTag << ")\n";
        return -1;
    }
    Contact c;
    c.tag = tag; c.masterSurfTag = masterSurfTag; c.slaveSurfTag = slaveSurfTag;
    c.kn = kn; c.kt = kt; c.mu = mu;
    theContacts.push_back(c);
    return 0;
}

int
LadrunoContactDomain::commit(void)
{
    // P1b: no path-dependent state yet (zero force). The counter lets the test
    // confirm the Domain::commit() hook actually drives us. P2 commits pair state.
    numCommits++;
    return 0;
}

int
LadrunoContactDomain::revertToLastCommit(void)
{
    numReverts++;
    return 0;
}
