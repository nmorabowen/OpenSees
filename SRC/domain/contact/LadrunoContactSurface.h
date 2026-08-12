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

// ADR-39 P1b — LadrunoContactSurface: a meshed contact surface, either a SLAVE
// node-set or a MASTER element-face / segment set. P1b stores only the topology
// (tags + segment connectivity); closest-point projection lands in P2.
//
// Lives in OPS_Domain (Domain-side state) alongside LadrunoContactDomain.
// See Ladruno_implementation/39_ladruno_contact_domain_adr.md + _adr39_p1b_design.md.

#ifndef LadrunoContactSurface_h
#define LadrunoContactSurface_h

#include <ID.h>

class Domain;

class LadrunoContactSurface
{
  public:
    enum Kind { SLAVE_NODES = 0, MASTER_SEGMENTS = 1, SLAVE_SEGMENTS = 2 };

    // SLAVE_NODES: nodeTags are the slave nodes.
    // MASTER_SEGMENTS: nodeTags is the flat connectivity, nodesPerSeg the stride
    //   (3 = tri, 4 = quad); segment s spans nodeTags[s*nodesPerSeg .. +stride).
    // SLAVE_SEGMENTS (ADR-41 C2.0): a FACETED slave surface — same flat-connectivity
    //   + nodesPerSeg layout as MASTER_SEGMENTS, but on the slave side. The mortar lane
    //   integrates the slave-slave matrix D over these facets (the node-set SLAVE_NODES
    //   carries no facet, so it cannot give D); the NTS lane keeps using SLAVE_NODES.
    LadrunoContactSurface(int tag, Kind kind, const ID &nodeTags, int nodesPerSeg = 0);
    ~LadrunoContactSurface();

    int getTag(void) const { return myTag; }
    Kind getKind(void) const { return myKind; }
    const ID &getNodeTags(void) const { return theNodeTags; }
    int getNumSegments(void) const;       // MASTER: count; SLAVE: numNodes
    int getNodesPerSeg(void) const { return nodesPerSeg; }

    // ADR-78 removal lane — drop nodes the ANALYSIS removed at runtime
    // (`recorder Collapse` / `remove node`), returning the number of node-tag
    // entries dropped. Nothing else in the engine may shrink a surface.
    //
    // This is only sound because a missing node can no longer mean a deck typo:
    // OPS_LadrunoContactSurface() now rejects a tag that does not exist at
    // DECLARATION time, so absence later can only be a runtime removal. Before
    // that check existed the two cases were indistinguishable at handle() time,
    // which is why ADR-78 P1 had to abort on both.
    //
    // Faceted kinds lose the WHOLE segment if any of its nodes is gone -- a
    // partial facet is not a facet, and the handler derives nSeg = size/nps, so
    // leaving a short tail would silently mis-stride every later segment.
    int pruneMissingNodes(Domain *dom);

  private:
    int myTag;
    Kind myKind;
    ID theNodeTags;
    int nodesPerSeg;
};

#endif
