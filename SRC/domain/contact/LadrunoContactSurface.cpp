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

// ADR-39 P1b — LadrunoContactSurface implementation. See header.

#include "LadrunoContactSurface.h"
#include <Domain.h>
#include <Node.h>

LadrunoContactSurface::LadrunoContactSurface(int tag, Kind kind,
                                             const ID &nodeTags, int npseg)
  : myTag(tag), myKind(kind), theNodeTags(nodeTags), nodesPerSeg(npseg)
{
}

LadrunoContactSurface::~LadrunoContactSurface()
{
}

int
LadrunoContactSurface::getNumSegments(void) const
{
    if ((myKind == MASTER_SEGMENTS || myKind == SLAVE_SEGMENTS) && nodesPerSeg > 0)
        return theNodeTags.Size() / nodesPerSeg;   // faceted: count facets
    return theNodeTags.Size();   // SLAVE_NODES: one entry per slave node
}

int
LadrunoContactSurface::pruneMissingNodes(Domain *dom)
{
    if (dom == 0)
        return 0;

    const int n = theNodeTags.Size();
    if (n == 0)
        return 0;

    const bool faceted =
        (myKind == MASTER_SEGMENTS || myKind == SLAVE_SEGMENTS) && nodesPerSeg > 0;
    const int stride = faceted ? nodesPerSeg : 1;

    // Whole strides only: for a faceted surface the stride is the facet, so one
    // missing corner retires the facet. For SLAVE_NODES the stride is 1, i.e. the
    // node itself.
    ID kept(0, n);
    int nkept = 0;
    for (int base = 0; base + stride <= n; base += stride) {
        bool alive = true;
        for (int k = 0; k < stride; k++) {
            if (dom->getNode(theNodeTags(base + k)) == 0) {
                alive = false;
                break;
            }
        }
        if (!alive)
            continue;
        for (int k = 0; k < stride; k++)
            kept[nkept++] = theNodeTags(base + k);
    }

    const int dropped = n - nkept;
    if (dropped > 0) {
        ID tight(nkept);
        for (int i = 0; i < nkept; i++)
            tight(i) = kept(i);
        theNodeTags = tight;
    }
    return dropped;
}
