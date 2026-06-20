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

// ADR-30 P1: LadrunoConstraintProjector — the M-orthogonal acceleration projector.
//
// For each connected-component constraint group it stores the stacked relation
// L (rows = all group DOFs, retained-first; cols = retained free DOFs), with the
// constrained rows = C. The enforced acceleration is the M-orthogonal projection
//     a_proj = L (Lt M L)^-1 Lt M a_raw           (M = diag of lumped masses)
// computed per group (tiny dense solve). Momentum conserved exactly (Lt M (a_proj
// - a_raw) = 0); equalDOF reduces to the mass-weighted average (Theory Eq. 28.1).
//
// Mass source (Gate-A R2): the per-equation lumped mass is read from the assembled
// DiagonalSOE diagonal — the exact M the integrator inverts — NOT from Node masses
// (which omit element-contributed mass). buildMass() must be called between the
// integrator's M-only formTangent() and its solve() (the DiagonalDirectSolver
// overwrites the diagonal with 1/aii on the factor pass).
//
// See Ladruno_implementation/30_..._adr.md and _adr30_p1_design.md (R1-R10).

#ifndef LadrunoConstraintProjector_h
#define LadrunoConstraintProjector_h

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>

class LinearSOE;

class LadrunoConstraintProjector
{
  public:
    LadrunoConstraintProjector();
    ~LadrunoConstraintProjector();

    // Called once by the handler (doneNumberingDOF) per group, with equation-index
    // data. Rows of L are ordered retained-first (nRet identity rows) then
    // constrained (the C rows). allEqn = equation number of each row (all >= 0).
    // delta = per-constrained-row offset (Uc0 - C*Ur0) for the IC check. fixedEqn =
    // slave equations whose every retained master is SP-fixed -> force-set to 0.
    void addGroup(const ID &allEqn, int nRet, const Matrix &L,
                  const Vector &delta, const ID &fixedEqn);

    int numGroups(void) const { return (int)groups.size(); }

    // Read diag(M) from the assembled Diagonal SOE and factor Lt M L per group.
    // Returns -1 (named error) on a non-Diagonal SOE, a massless retained
    // direction, or a singular Lt M L. Sets the mass-ready flag.
    int buildMass(LinearSOE *theSOE);

    // In-place M-orthogonal projection of the acceleration vector a (size = numEqn).
    void project(Vector &a);

    // IC compliance of committed (u,v) against every group, in the MP increment
    // convention (u_c - C u_r = delta; v_c - C v_r = 0). Returns the number of
    // violated constrained rows (0 = compliant); prints a named diagnostic per
    // violation. Never aborts itself — the caller decides.
    int checkIC(const Vector &u, const Vector &v, double icTol) const;

    // -projectICs: snap the constrained components of u and v to the master-implied
    // compliant values (u_c = C u_r + delta; v_c = C v_r; orphaned slaves -> 0).
    void snapICs(Vector &u, Vector &v) const;

    // IC policy set by the handler. enforceIC() is the single call the integrator
    // makes each domainChanged(): if -projectICs it snaps (u,v) compliant and returns
    // 0; otherwise it checks and returns the violation count (>0 -> caller aborts).
    void setICPolicy(bool projectICs_, double icTol_) { projectICs = projectICs_; icTol = icTol_; }
    int enforceIC(Vector &u, Vector &v) const;

    bool isMassReady(void) const { return massReady; }
    void invalidateMass(void) { massReady = false; }

    // Tie-force query (ADR-30 P3): the constraint force f_tie = M (a_raw - a_proj) per
    // DOF, cached during the last project(). Indexed by global equation number; 0 on a
    // DOF not in any group. Returns 0 if no projection has run yet.
    double tieForceAtEqn(int eqn) const;
    const Vector &getTieForces(void) const { return tieForce; }

  private:
    struct Group {
        ID allEqn;      // eqn of each L row (retained-first), all >= 0
        int nRet;       // number of retained rows/cols (L top block = I_nRet)
        Matrix L;       // (nRet+nCon) x nRet
        Vector delta;   // nCon offsets for the constrained rows (IC check)
        Vector m;       // per-row lumped mass (filled by buildMass)
        Matrix LtML;    // nRet x nRet, factored-on-demand copy used in project()
        ID fixedEqn;    // slave eqns force-set to 0 (orphaned by SP-fixed master)
        Group() : nRet(0) {}
    };
    std::vector<Group> groups;
    bool massReady;
    bool projectICs;
    double icTol;
    Vector tieForce;   // P3: cached f_tie = M(a_raw - a_proj), sized to numEqn at project()
};

#endif
