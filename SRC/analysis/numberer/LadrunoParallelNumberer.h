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
** ****************************************************************** */

// Ladruno — LadrunoParallelNumberer (ADR-74)
//
// File: SRC/analysis/numberer/LadrunoParallelNumberer.h
//
// Written: Ladruno / nmora
// Created: 2026-07
//
// Description: The O(V) replacement engine for ParallelNumberer's O(V²)
//   gather/merge (the 16.35 h first-step stall on 19.18 M-node np240 runs —
//   measured dc.numberDOF ~ N^2.01, ADR-74 G0). Fork-native subclass in the
//   ≥33000 classTag band via the ledgered ParallelNumberer.h promotion;
//   upstream classes stay behavior-untouched. Deck verbs:
//     numberer LadrunoParallelRCM     (RCM on the merged graph — G1
//                                      bit-identical to stock ParallelRCM)
//     numberer LadrunoParallelPlain   (subdomain ordering, no RCM pass —
//                                      the recommended explicit-lane verb;
//                                      G1b bijection + same-physics)
//
//   Phase status (ADR-74 §Implementation plan):
//     N1 ✓: delegate mode, proven bit-identical to stock (13/13, PR #591).
//     N2 (this, T0): own numberDOF(int) — O(1) dedup in both merge passes
//        (RefIndex dense/hash), per-subgraph tag map, seen-flag plain branch
//        (fixes the upstream tag-0 quirk => G1b not bit-identity), ref<0
//        Lagrange guard, constrainedNode->MP index for the -4 fixup,
//        dc.n.gather/merge/order/scatter profiler sub-scopes. RCM verb stays
//        bit-identical to stock (G1). numberDOF(ID&) hard-errors: the
//        upstream variant is broken (ParallelNumberer.cpp:507+, numbering
//        call commented out) and is SP-lane-only, which ADR-74 scopes out.
//     N3 (T1): owned dense tag->Vertex* + Vertex::addEdge direct + dense-
//        storage merged graph (kills the std::map residual).
//
// Theory + gates: Ladruno_implementation/74_ladruno_parallel_numberer_adr.md

#ifndef LadrunoParallelNumberer_h
#define LadrunoParallelNumberer_h

#include <ParallelNumberer.h>

class GraphNumberer;

class LadrunoParallelNumberer : public ParallelNumberer
{
  public:
    // RCM (or any GraphNumberer) form — the LadrunoParallelRCM verb.
    // Ownership matches the base: theGraphNumberer is deleted by the base dtor.
    LadrunoParallelNumberer(GraphNumberer &theGraphNumberer);

    // no-numberer form (subdomain ordering) — the LadrunoParallelPlain verb
    LadrunoParallelNumberer();

    ~LadrunoParallelNumberer();

    // T0: the O(V) gather/merge/order/scatter (see the .cpp header).
    int numberDOF(int lastDOF = -1);

    // Upstream numberDOF(ID&) is broken (see header comment) — hard error
    // instead of faithfully wrapping code that never worked.
    int numberDOF(ID &lastDOFs);
};

#endif
