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
// References:
//   Wriggers, P. (2006). Computational Contact Mechanics, 2nd ed. Springer,
//     Berlin. doi:10.1007/978-3-540-32609-0
//   Laursen, T.A. (2002). Computational Contact and Impact Mechanics.
//     Springer, Berlin. doi:10.1007/978-3-662-04864-1
//   Hallquist, J.O. (2006). LS-DYNA Theory Manual, Livermore Software
//     Technology Corporation, Section 26 (Contact-Impact Algorithm).
//
// Description: LadrunoContactDomain, the contact engine attached to a Domain.
// It owns the contact surfaces and the contact definitions, and the
// path-dependent pair state, so that the assembly-side adapters
// (LadrunoContactFE, rebuilt in every ConstraintHandler::handle()) are
// stateless views.
//
// Lifetime: a Domain owns at most one engine (none means no cost and stock
// behaviour). The engine survives domainChanged() (which clears the
// AnalysisModel, not the Domain) and is deleted by Domain::clearAll() (wipe)
// and by ~Domain. Domain::commit(), revertToLastCommit() and revertToStart()
// call the matching engine methods, so path-dependent contact state commits
// and reverts with the rest of the model whatever the integrator.

#ifndef LadrunoContactDomain_h
#define LadrunoContactDomain_h

#include <vector>
#include <map>
#include <set>

class LadrunoContactSurface;
class Domain;

class LadrunoContactDomain
{
  public:
    LadrunoContactDomain();
    ~LadrunoContactDomain();

    // --- model definition (contactSurface / contact commands) ---
    int addSurface(LadrunoContactSurface *surf);          // takes ownership
    LadrunoContactSurface *getSurface(int tag) const;

    // Reconcile the declared contact with a domain the analysis has been
    // shrinking (`remove node`, element removal). Called at the start of every
    // handle(): removed nodes are pruned from every surface, and an interaction
    // whose master or slave surface is left empty is retired; both are reported
    // on opserr. Returns the number of interactions retired by this call.
    int pruneRemovedNodes(Domain *dom);

    // A node-to-segment contact between a MASTER_SEGMENTS surface and a
    // SLAVE_NODES surface; the handler builds one LadrunoContactFE per
    // (slave node, candidate master segment) pair.
    //   outward:  optional direction toward the slave's allowed half-space,
    //             used to orient the derived segment normal; null = derived
    //             per pair from the slave and segment reference positions.
    //   knAuto:   kn is sized per pair at handle() time from the stiffness of
    //             the solid element that owns the master segment (kn is then
    //             a placeholder).
    //   cellFrac: broad-phase cell size as a multiple of the median segment
    //             diagonal (1.0 default; a huge value gives one bucket, i.e.
    //             brute force).
    //   consistentNormal: add the geometric part kn gN d^2gN/du^2 of the
    //             consistent normal tangent (curved or large-sliding masters).
    int addContact(int tag, int masterSurfTag, int slaveSurfTag,
                   double kn, const double *outward = 0,
                   bool knAuto = false, double cellFrac = 1.0,
                   bool consistentNormal = false);
    int getNumContacts(void) const { return (int)theContacts.size(); }

    struct Contact {
        // Set when the analysis removed every node of one of the surfaces of
        // this interaction. A retired interaction is skipped by handle() with a
        // named message; it cannot hide a typo, because a node tag that does
        // not exist is refused when the surface is declared.
        bool retired = false;
        int tag, masterSurfTag, slaveSurfTag;
        double kn;
        bool knAuto;            // kn sized from the master solid element
        bool hasOutward;        // an explicit orientation direction was given
        double outward[3];      // orientation direction toward the allowed half-space
        double cellFrac;        // bucket-sort cell = cellFrac * median segment diagonal
        bool consistentNormal;  // assemble the geometric normal tangent
    };
    const Contact &getContact(int i) const { return theContacts[i]; }

    // --- node-to-segment normal force readout (ladrunoContactForce) ---
    // Each segment adapter writes its normal penalty force tn = kn <-gap>+ into a
    // (contactTag, slaveTag, segIndex) slot at every residual evaluation, and
    // writes 0 when the pair is not active, so after convergence the slots hold
    // the converged forces of the pairs that are in contact. getNtsForce sums the
    // slots of one slave node. The slots are cleared at every handle(). No effect
    // on the assembled residual or tangent.
    void   clearNtsForce(void);
    void   setNtsForce(int contactTag, int slaveTag, int segIndex, double tn);
    double getNtsForce(int slaveTag) const;   // sum over (contactTag, slaveTag, *)

    // --- lifecycle (driven by Domain::commit / revertToLastCommit / revertToStart) ---
    // augmenting: the Domain's within-step augmentation flag, passed for
    // contact formulations that update multipliers on commit.
    int commit(bool augmenting = false);
    int revertToLastCommit(void);
    // Drop all path-dependent state so that Domain::revertToStart() returns
    // contact to the state of a fresh run; the definitions are kept.
    int revertToStart(void);
    int getNumCommits(void) const { return numCommits; }
    int getNumReverts(void) const { return numReverts; }

  private:
    std::vector<LadrunoContactSurface *> theSurfaces;   // owned
    std::vector<Contact> theContacts;
    // The contact tag is the leading key of every path-state store, so tags
    // must be unique across all contact definitions (checked by every add*).
    bool contactTagInUse(int tag) const;
    int numCommits;
    int numReverts;

    // (contactTag, slaveTag, segIndex) key of the pair state stores.
    struct PairKey {
        int c, s, g;
        bool operator<(const PairKey &o) const {
            if (c != o.c) return c < o.c;
            if (s != o.s) return s < o.s;
            return g < o.g;
        }
    };
    std::map<PairKey, double> theNtsForce;              // per-pair normal force readout
};

#endif
