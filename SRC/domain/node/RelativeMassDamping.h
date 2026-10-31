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

// Ladruno (relativeMassDamping): mass-proportional Rayleigh damping relative to a
// prescribed reference velocity (the free field of an absolute-coordinate SSI model).
//
// With a reference attached, a node's mass-proportional damping force is
//
//     f_D = alphaM * M * (v - v_ref(t)),     v_ref(t)_i = series_i(t)  (0 if none)
//
// instead of alphaM * M * v. v_ref is a prescribed function of the domain time
// only, so the damping matrix (Node::getDamp = alphaM * M) and therefore every
// tangent are unchanged; only the residual, the Rayleigh / reaction responses
// and the energy bookkeeping see the reference term.
//
// Objects are owned by a process-wide registry (one per tag), deleted by
// OPS_clearAllRelativeMassDampings() from Domain::clearAll() (wipe). Nodes hold
// a non-owning pointer. Absent any reference, every code path is the original.
//
// Command (classic Tcl, interpreter Tcl, Python; the name is OPS_RelativeMassDampingCmdName):
//   relativeMassDamping <tag> -vel <dof> <tsTag> [-vel <dof> <tsTag> ...]
//                       ( -region <regTag> ... | -node <nodeTag> ... | -allMassNodes )
// Attaches the reference to the listed nodes that carry a nodal mass on this
// process. Returns { nAttached nAlphaZero sum(alphaM*m_dof1) sum(alphaM*m_dof2) ... }
// for the referenced dofs, in -vel order, for the audit.

#ifndef RelativeMassDamping_h
#define RelativeMassDamping_h

#include <Vector.h>

class TimeSeries;

class RelativeMassDamping
{
  public:
    // series[i] is the series of dof dofs(i) (0-based); copies are taken.
    RelativeMassDamping(int tag, int numRef, const int *dofs, TimeSeries **series);
    ~RelativeMassDamping();

    int getTag(void) const {return tag;}

    // v_ref at domain time t, sized numDOF (dofs beyond numDOF ignored)
    const Vector &getVelocity(double time, int numDOF);

    // the time of the last evaluation (for diagnostics)
    double getLastTime(void) const {return lastTime;}

  private:
    int tag;
    int numRef;
    int *dofs;
    TimeSeries **series;
    double *factors;      // cached series values at lastTime
    double lastTime;
    bool   haveCache;
    Vector result;        // scratch, resized on demand
};

// registry
int  OPS_addRelativeMassDamping(RelativeMassDamping *theRef);
RelativeMassDamping *OPS_getRelativeMassDamping(int tag);
int  OPS_numRelativeMassDampings(void);
void OPS_clearAllRelativeMassDampings(void);

// the command name: the ONE place to rename it (all registrations and messages read it)
extern const char *OPS_RelativeMassDampingCmdName;

// command entry point (elementAPI cursor); returns 0 on success, -1 on error
int OPS_RelativeMassDamping(void);

#endif
