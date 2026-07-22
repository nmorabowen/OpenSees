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

// Ladruno — LadrunoParallelNumberer (ADR-74). See the header for the phase
// plan; this file is the N1 DELEGATE — provably bit-identical to stock before
// T0/T1 change the algorithm underneath.

#include <LadrunoParallelNumberer.h>
#include <classTags.h>
#include <OPS_Globals.h>

LadrunoParallelNumberer::LadrunoParallelNumberer(GraphNumberer &theGraphNumberer)
  : ParallelNumberer(NUMBERER_TAG_LadrunoParallelNumberer, theGraphNumberer)
{

}

LadrunoParallelNumberer::LadrunoParallelNumberer()
  : ParallelNumberer(NUMBERER_TAG_LadrunoParallelNumberer)
{

}

LadrunoParallelNumberer::~LadrunoParallelNumberer()
{
  // base dtor owns theNumberer + the channel array
}

int
LadrunoParallelNumberer::numberDOF(int lastDOF)
{
  // N1 delegate: the G1 null-test of the oracle AND of this class's wiring.
  // T0 (N2) replaces this with the O(V) merge; the gate stays bit-identity.
  return this->ParallelNumberer::numberDOF(lastDOF);
}

int
LadrunoParallelNumberer::numberDOF(ID &lastDOFs)
{
  // Upstream ParallelNumberer::numberDOF(ID&) is non-functional (the numbering
  // call is commented out and the recv layout mismatches — see ADR-74 §Where).
  // It is reachable only from the SP/DomainDecomposition lane, which ADR-74
  // scopes out. Fail loudly rather than wrap code that never worked.
  opserr << "ERROR LadrunoParallelNumberer::numberDOF(ID&) — the lastDOFs "
         << "variant is unsupported (upstream implementation is broken; "
         << "SP lane out of scope, ADR-74). Use numberDOF(int).\n";
  return -1;
}
