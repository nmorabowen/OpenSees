/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Developed by: N. Mora Bowen, P. Palacios, J.A. Abell (Ladruño project)
//
// Description: LadrunoContactDomain, see LadrunoContactDomain.h.

#include "LadrunoContactDomain.h"
#include "LadrunoContactSurface.h"
#include <OPS_Globals.h>
#include <cmath>
#include <algorithm>
#include <Domain.h>
#include <Node.h>
#include <ID.h>

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
    if (getSurface(surf->getTag()) != 0) {
        opserr << "WARNING LadrunoContactDomain::addSurface() - surface tag "
               << surf->getTag() << " already exists\n";
        return -1;
    }
    // The handler derives nSeg = size/nodesPerSeg, so a faceted surface must
    // hold a whole number of segments; a trailing partial segment would be
    // dropped without a message.
    if (surf->getKind() != LadrunoContactSurface::SLAVE_NODES) {
        int nps = surf->getNodesPerSeg();
        int n   = surf->getNodeTags().Size();
        if (nps < 3 || nps > 4 || (n % nps) != 0) {
            opserr << "WARNING LadrunoContactDomain::addSurface() - surface tag "
                   << surf->getTag() << ": " << n << " connectivity entries is not a whole "
                      "number of " << nps << "-node segments (nodesPerSeg must be 3 or 4)\n";
            return -1;
        }
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

bool
LadrunoContactDomain::contactTagInUse(int tag) const
{
    for (size_t i = 0; i < theContacts.size(); i++)
        if (theContacts[i].tag == tag) return true;
    return false;
}

int
LadrunoContactDomain::addContact(int tag, int masterSurfTag, int slaveSurfTag,
                                 double kn, const double *outward,
                                 bool knAuto, double cellFrac, bool consistentNormal)
{
    if (getSurface(masterSurfTag) == 0 || getSurface(slaveSurfTag) == 0) {
        opserr << "WARNING LadrunoContactDomain::addContact() - master/slave surface "
                  "not defined (master " << masterSurfTag << ", slave " << slaveSurfTag << ")\n";
        return -1;
    }
    if (contactTagInUse(tag)) {
        // the tag keys the pair state stores; a duplicate would alias two interactions
        opserr << "WARNING LadrunoContactDomain::addContact() - contact tag " << tag
               << " already exists (contact tags key the pair state)\n";
        return -1;
    }
    if (kn < 0.0) {
        // kn < 0 is an attractive contact with a negative-definite tangent.
        // kn == 0 is accepted here (a topology-only definition); handle() refuses
        // a segment contact without a positive or automatic penalty.
        opserr << "WARNING LadrunoContactDomain::addContact() - penalty kn must be >= 0 (got "
               << kn << ")\n";
        return -1;
    }
    Contact c;
    c.tag = tag; c.masterSurfTag = masterSurfTag; c.slaveSurfTag = slaveSurfTag;
    c.kn = kn;
    c.knAuto = knAuto;
    c.hasOutward = (outward != 0);
    for (int d = 0; d < 3; d++) c.outward[d] = (outward != 0) ? outward[d] : 0.0;
    c.cellFrac = (cellFrac > 0.0) ? cellFrac : 1.0;
    c.consistentNormal = consistentNormal;
    theContacts.push_back(c);
    return 0;
}

void
LadrunoContactDomain::clearNtsForce(void)
{
    theNtsForce.clear();
}

void
LadrunoContactDomain::setNtsForce(int contactTag, int slaveTag, int segIndex, double tn)
{
    PairKey k; k.c = contactTag; k.s = slaveTag; k.g = segIndex;
    theNtsForce[k] = tn;        // overwritten at every residual evaluation
}

double
LadrunoContactDomain::getNtsForce(int slaveTag) const
{
    double f = 0.0;
    for (std::map<PairKey, double>::const_iterator it = theNtsForce.begin();
         it != theNtsForce.end(); ++it)
        if (it->first.s == slaveTag) f += it->second;
    return f;
}

int
LadrunoContactDomain::commit(bool)
{
    numCommits++;
    return 0;
}

int
LadrunoContactDomain::revertToLastCommit(void)
{
    numReverts++;
    return 0;
}

int
LadrunoContactDomain::revertToStart(void)
{
    // Every slot is recreated zeroed at the next handle() or residual
    // evaluation, i.e. in the state of a fresh run.
    theNtsForce.clear();
    return 0;
}

int
LadrunoContactDomain::pruneRemovedNodes(Domain *dom)
{
    if (dom == 0)
        return 0;

    // 1. drop the nodes the analysis removed, keeping facets whole
    for (size_t i = 0; i < theSurfaces.size(); i++) {
        LadrunoContactSurface *s = theSurfaces[i];
        if (s == 0)
            continue;
        int before = s->getNodeTags().Size();
        int dropped = s->pruneMissingNodes(dom);
        if (dropped > 0) {
            opserr << "LadrunoContactDomain: contactSurface " << s->getTag()
                   << " lost " << dropped << " of " << before
                   << " node entries removed by the analysis; "
                   << s->getNodeTags().Size() << " remain\n";
        }
    }

    // 2. retire every interaction left without a usable surface
    int retiredNow = 0;
    for (size_t i = 0; i < theContacts.size(); i++) {
        Contact &c = theContacts[i];
        if (c.retired)
            continue;
        LadrunoContactSurface *ms = getSurface(c.masterSurfTag);
        LadrunoContactSurface *ss = getSurface(c.slaveSurfTag);
        if (ms == 0 || ss == 0 ||
            ms->getNodeTags().Size() == 0 || ss->getNodeTags().Size() == 0) {
            c.retired = true;
            retiredNow++;
            opserr << "LadrunoContactDomain: RETIRING contact " << c.tag
                   << " -- the analysis removed every node of its master or slave "
                      "surface, so the interaction no longer has a body to act on. "
                      "The run continues without it.\n";
        }
    }
    return retiredNow;
}
