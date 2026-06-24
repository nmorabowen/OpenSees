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

// ADR-39 — LadrunoContactFE: an FE_Element adapter that injects contact
// contributions into the assembly without being backed by a Domain Element.
//
// The contact narrow phase is self-contained here: getResidual()/getTangent()
// read node trial displacements (like a real Element) and assemble F_c / K_c.
// Explicit integrators consume only getResidual() (formEleTangent is mass-only
// under CentralDifferenceLadruno); implicit consume both (tangent scaled by the
// integrator's c1). The adapter is a STATELESS VIEW — all path-dependent pair
// state lives on the Domain-owned LadrunoContactDomain (P1b), so the adapter can
// be destroyed/rebuilt every handle() with no state loss.
//
// Three modes, shipped incrementally:
//   EMPTY       (P1a, #345) — size-0 getID/getDOFtags + ZERO force; legacy/null path proving the
//                             injection plumbing is graph-neutral (bitwise-identical to no-contact).
//   RIGID_PLANE (P2a, #346) — slave node vs a rigid analytical plane (penalty normal).
//   SEGMENT     (P2b/P3/P3.5, #354/#360/#361) — faceted NTS penalty + Coulomb friction + consistent
//                             tangent, real connectivity FE_Element(tag, 1+nps, 3*(1+nps)).
// With myEle==0 (the EMPTY mode) the base getResidual/getTangent exit(-1) and the base zero/add
// helpers early-return, so this subtype MUST own its own buffers and override both.
//
// See Ladruno_implementation/39_ladruno_contact_domain_adr.md + _adr39_p1_design.md.

#ifndef LadrunoContactFE_h
#define LadrunoContactFE_h

#include <FE_Element.h>
#include <Vector.h>
#include <Matrix.h>

class Integrator;
class Node;
class Domain;

class LadrunoContactFE : public FE_Element
{
  public:
    // P1a: empty-connectivity adapter (numDOF_Group = 0, ndof = 0).
    LadrunoContactFE(int tag);
    // P2a: rigid analytical plane vs ONE slave node. Connectivity = {slave} (the
    // first ndm translational DOFs). p0 = a point on the plane, n = outward unit
    // normal (toward the slave's allowed half-space), kn = penalty stiffness.
    // D2: muc = viscous normal-stabilization coefficient (p_visc = muc*gap_rate); 0 ⇒ no viscous term.
    LadrunoContactFE(int tag, Node *slaveNode, int ndm,
                     const double p0[3], const double n[3], double kn, double muc = 0.0);
    // P2b: faceted node-to-segment penalty contact vs ONE master segment (tri-3 or
    // quad-4, nps nodes). Connectivity = {slave} ∪ {segment nodes}; the outward
    // normal is DERIVED per evaluation + oriented toward orientDir (the slave's
    // allowed half-space — winding-immune, design-gate BLOCKER-1). 3D only (ndm==3).
    //
    // P3 friction: kt (tangential penalty), mu (Coulomb). When mu>0 the adapter adds
    // the Coulomb friction traction, keeping its path-dependent plastic-slip state on
    // the Domain-owned LadrunoContactDomain (looked up lazily via theDomain each
    // getResidual — wipe deletes the engine, so the adapter must NOT cache its ptr),
    // keyed (contactTag, slaveTag, segIndex). mu<=0 ⇒ byte-identical to frictionless
    // P2b (no slot touch). kt/mu/engine args default to the frictionless path.
    // P3.5 implicit friction tangent: consistentTan=false (default) assembles the
    // SYMMETRIC tangent (drop the d_TN⊗n column ⇒ solver-safe on any system, superlinear);
    // true assembles the full non-symmetric consistent tangent (quadratic, needs
    // FullGeneral/UmfPack). Only ever assembled under IMPLICIT integrators (addKtToTang);
    // explicit CDL routes to addMtoTang (no-op) ⇒ P3 explicit byte-identical.
    // B3 (P2b-2c) consistentNormal: assemble the consistent ∂n/∂u geometric NORMAL
    // tangent K_geom = kn·gN·∂²gN/∂u² (the curvature / large-sliding block deferred since
    // P2b). false (default) ⇒ the shipped kn·BᵀB main term only (byte-identical, the flat/
    // fixed-master EXACT tangent). true ⇒ + the geometric block ⇒ quadratic Newton on
    // CURVED / large-sliding interfaces. SYMMETRIC (Hessian of the scalar gap) ⇒ solver-safe.
    LadrunoContactFE(int tag, Node *slaveNode, Node **segNodes, int nps, double kn,
                     const double orientDir[3], double kt = 0.0, double mu = 0.0,
                     Domain *theDomain = 0, int contactTag = 0, int segIndex = 0,
                     bool consistentTan = false, double muc = 0.0,
                     bool consistentNormal = false);
    // C2.1/C2.2 (ADR-41): clipped-GP MORTAR contact, ONE slave facet vs ONE master facet
    // (tri-3 / quad-4). Connectivity = {slave facet nodes} ∪ {master facet nodes},
    // FE_Element(tag, nps_s+nps_m, 3*(nps_s+nps_m)). getResidual integrates the pair via
    // LadrunoMortarKernel::integratePair (D,M,g̃ at the trial config) and assembles the
    // augmented force F^s=−(D·p)n, F^m=+(Mᵀ·p)n with the per-GLOBAL-slave-node pressure
    // p_I = min(0, λ_I + epsN·ḡ_I). C2.2: λ_I and the GLOBAL weighted gap ḡ_I live on the
    // Domain-owned LadrunoContactDomain (looked up lazily via theDomain each getResidual —
    // wipe deletes the engine, so the adapter must NOT cache its ptr), keyed (contactTag,
    // slaveNodeTag); the adapter reports its per-facet g̃_I^facet/a_I^facet into that node slot
    // (accumulateMortarGap) then reads the running global gap. theDomain==0 ⇒ the C2.1 penalty
    // fallback (λ≡0, facet-local gap). addKtToTang assembles K_c=epsN·B̃ᵀdiag(act/a_global)B̃⊗(n⊗n).
    // C3.1 friction: mu/epsT(in the kt slot)/cohesion/tauMax add the Coulomb/Tresca tangential
    // traction (shipped LadrunoFrictionKernel). All ≤ 0 ⇒ byte-identical to the frictionless C2
    // path (no slot touch — the NTS P3 `mu>0` short-circuit, generalized to the unified cone).
    // consistentTan = the non-symmetric Coulomb friction tangent (C3.2; false = symmetric).
    // C4 isTie: a MESH-TIE pair (a permanent bond — the zero-gap limit of contact). When true the
    // adapter skips the normal KKT + friction blocks entirely and assembles the FULL 3-vector tie
    // force t_I = λ_tie,I + epsTie·(r_I/a_I) (no clamp) via D/−M + the SPD tangent epsTie·B̃ᵀB̃⊗I₃.
    // epsN carries epsTie; mu/cohesion/tauMax are refused with -tie upstream. theDomain==0 ⇒ a pure
    // penalty tie (λ_tie≡0). Mutually exclusive with the friction args (a tie has no cone).
    LadrunoContactFE(int tag, Node **slaveNodes, int nps_s, Node **masterNodes, int nps_m,
                     double epsN, const double orientDir[3], int contactTag = 0,
                     int slaveFacetIndex = 0, Domain *theDomain = 0,
                     double mu = 0.0, double epsT = 0.0, double cohesion = 0.0,
                     double tauMax = 0.0, bool consistentTan = false, bool isTie = false,
                     double muc = 0.0);   // D2.2: viscous normal stabilization on the mortar contact
    ~LadrunoContactFE();

    // self-owned buffers (base buffers are unavailable when myEle == 0)
    const Vector &getResidual(Integrator *theIntegrator);
    const Matrix &getTangent(Integrator *theIntegrator);

    // getTangent routes through the integrator's formEleTangent so the INTEGRATOR
    // decides what to assemble (CDL -> addMtoTang only -> no contact stiffness in
    // the explicit mass matrix; Newmark -> addKtToTang(c1) -> c1*K_c; statics ->
    // addKtToTang(1) -> K_c). These overrides feed MY tang buffer (base helpers
    // early-return when myEle==0).
    void zeroTangent(void);
    void addKtToTang(double fact = 1.0);  // += fact * K_c when the pair is active
    void addKiToTang(double fact = 1.0);  // initial-stiffness path: K_initial == K_c
                                          // here (flat rigid plane, n constant), so
                                          // mirror addKtToTang -> Newton -initial /
                                          // ModifiedNewton -initial keep the contact
                                          // stiffness instead of silently dropping it.
    void addCtoTang(double fact = 1.0);   // contact has no damping -> no-op
    void addMtoTang(double fact = 1.0);   // contact pairs carry no mass -> no-op

    // With myEle==0 the base force-vector helpers print a WARNING and return the
    // shared error vector (FE_Element.cpp). They are reachable off the assembly
    // hot path (e.g. modal damping doMv -> getM_Force when mass is non-diagonal),
    // so override them to a clean size-0 return. (P2: real K/M/C force variants.)
    const Vector &getTangForce(const Vector &x, double fact = 1.0);
    const Vector &getK_Force(const Vector &x, double fact = 1.0);
    const Vector &getKi_Force(const Vector &x, double fact = 1.0);
    const Vector &getC_Force(const Vector &x, double fact = 1.0);
    const Vector &getM_Force(const Vector &x, double fact = 1.0);

  private:
    // gap n.(x_s - p0) at the slave's current trial position; <0 = penetration
    double rigidPlaneGap(void) const;

    // P2b: project the slave onto the master segment at the current trial config,
    // derive the oriented normal + gap + shape weights. Returns the ACTIVE flag
    // (true = penetrating + in-bounds). Fills gap(<0), n[3], N[nps], and the
    // assembled gap operator B[ndof] = [ nᵀ | −N_i nᵀ ] over [slave, seg nodes].
    // P3: if gTvec != 0, also fills the tangential relative-position vector
    // (x_s − x̄)_tangential (slip measure for the friction return map).
    bool segmentActive(double &gap, double n[3], double N[4], double *B,
                       double *gTvec = 0) const;

    // P3.5: scatter the 3×3 friction tangent block K_ss (LadrunoContactKernel::
    // frictionTangentBlock) into `tang` via G = [I | −N_i I]: block (a,b) = w_a w_b K_ss
    // with w = [1, −N_0, …, −N_nps]. consistent ⇒ include the non-symmetric d_TN⊗n.
    void addFrictionTang(double fact, const double n[3], const double N[4], double tn,
                         const double gTeff[3], const double gpT[3], bool consistent);

    // B3 (P2b-2c): assemble the consistent ∂n/∂u geometric NORMAL tangent block
    // K_geom = kn·gN·∂²gN/∂u² into `tang` (SEGMENT mode). Re-projects the current config
    // (deterministic ⇒ the SAME ξ̄/n/gap the residual used) and calls the oracle-validated
    // LadrunoContactProjection::normalGeomTangent. Active mask = penetrating + in-bounds;
    // for a flat facet the slave block is identically 0 (byte-identity contract). SYMMETRIC.
    void addNormalGeomTang(double fact);

    // C2.1: integrate the mortar facet pair at the current trial config. Fills D,M,g̃
    // (LadrunoMortarKernel::integratePair) + the per-facet master normal n. Returns true
    // if the overlap is non-empty (status 0); false ⇒ no contribution this evaluation.
    bool mortarActive(double D[4][4], double M[4][4], double g[4], double n[3]) const;
    // C2.1: assemble the mortar penalty tangent K_c = epsN·B̃ᵀ diag(act/a) B̃ ⊗ (n⊗n)
    // into `tang` (material/penalty only — geometric ∂{D,M,n}/∂u deferred). Shared by
    // addKtToTang / addKiToTang (the penalty K_initial == K_current). Same active mask
    // as the residual. fact = the integrator's c1 (or 1 for statics).
    // initialStiff=true (addKiToTang path) forces the friction STICK tangent kt·P_t (SPD ⇒ a
    // Modified/Initial-Newton contraction), instead of the rank-deficient slip tangent — mirrors
    // the SEGMENT addKiToTang rule (gate Q5). Default false (addKtToTang = current slip/stick).
    void addMortarTang(double fact, bool initialStiff = false);

    Vector resid;   // size-0 in P1a; ndm in P2a; ndm*(1+nps) in P2b
    Matrix tang;    // 0x0 in P1a; ndm x ndm in P2a; ndof x ndof in P2b

    // P2a rigid-plane binding (mode = RIGID_PLANE); unused/zero in P1a
    enum Mode { EMPTY = 0, RIGID_PLANE = 1, SEGMENT = 2, MORTAR = 3 };
    Mode mode;
    Node *theSlave;
    int ndm;
    double planeP0[3];
    double planeN[3];
    double kn;

    // P2b SEGMENT binding (mode == SEGMENT)
    Node *segNode[4];   // master segment nodes (tri-3 → 3, quad-4 → 4)
    int nps;            // nodes per segment
    double orientDir[3];// fixed direction toward the slave's allowed half-space
                        // (the derived normal is flipped to satisfy n·orientDir>0)

    // P3 friction binding (active only in SEGMENT mode with mu>0)
    double kt;          // tangential penalty
    double mu;          // Coulomb friction coefficient (<=0 ⇒ frictionless P2b path)
    double muc;         // D2: viscous normal-stabilization coeff (p_visc = muc*gap_rate; 0 ⇒ off).
                        // RIGID_PLANE/SEGMENT (D2.1) + MORTAR contact (D2.2, NOT ties); force-only under
                        // CDL (explicit), force + a C=muc·(B/B̃)ᵀ(B/B̃) damping tangent (addCtoTang) under implicit.
    Domain *theDomain;  // for the LAZY engine re-fetch (wipe deletes the engine, so
                        // we must not cache LadrunoContactDomain*); null ⇒ no friction
    int contactTag;     // friction-state key: contact definition tag ...
    int segIndex;       // ... and the GLOBAL master-segment ordinal (rebuild-stable)
    bool consistentTan; // P3.5: include the non-symmetric d_TN⊗n column (false ⇒ symmetric)
    bool consistentNormal; // B3 (P2b-2c): SEGMENT — add the consistent ∂n/∂u geometric normal
                        // tangent (kn·gN·∂²gN/∂u²). false ⇒ shipped kn·BᵀB only (byte-identical).

    // C2.1 MORTAR binding (mode == MORTAR). The penalty epsN rides `kn`. The slave-facet
    // ordinal `slaveFacetIndex` (+ contactTag) keys the per-node λ_N state in C2.2.
    Node *mortarSlave[4];   // slave facet nodes  (tri-3 → 3, quad-4 → 4)
    Node *mortarMaster[4];  // master facet nodes
    int npsS, npsM;         // slave / master nodes-per-facet
    int slaveFacetIndex;    // GLOBAL slave-facet ordinal (rebuild-stable; C2.2 λ_N key)

    // C3.1 MORTAR friction (active in MORTAR mode when mu>0 ∨ cohesion>0 ∨ tauMax>0). epsT (the
    // tangential penalty) rides `kt`; mu reuses the friction `mu` member. cohesion/tauMax complete
    // the unified cone cap = min(μN+c, τmax). consistentTan reuses the friction member.
    double mortarCohesion;  // adhesive intercept c
    double mortarTauMax;    // Tresca shear cap (≤0 ⇒ no upper cap)
    bool   isTie;           // C4: MESH-TIE pair (full 3-vec r→0, no clamp/friction); kn carries epsTie

    // C3.1: assemble the mortar Coulomb/Tresca friction force into `resid` (called from the MORTAR
    // getResidual after the normal block). p_normal[I] = the per-node normal pressure (≤0; the
    // C2 value); for each in-contact node it builds the LOCAL weighted tangential gap, runs the
    // shipped return map (engagement origin + committed slip on the Domain), and scatters the
    // tangential traction via D/M exactly like the normal force. cd = the (non-null) engine.
    void addMortarFriction(const double D[4][4], const double M[4][4], const double n[3],
                           const double p_normal[4], class LadrunoContactDomain *cd);

    // C4: assemble the MESH-TIE force into `resid`. For each slave node I build the FULL 3-vec
    // weighted relative DISPLACEMENT r_I = Σ_J D_IJ u_s,J − Σ_K M_IK u_m,K (from getTrialDisp — the
    // bond exists from the as-built config, so no gT0), form the tie traction t_I = λ_tie,I +
    // epsTie·(r_I/a_I) (NO clamp), and scatter it via D/−M like the normal force (f^s_K = −Σ D_KI t_I,
    // f^m_L = +Σ M_IL t_I). Reports r_I into the Domain GLOBAL accumulator (for the commit Uzawa +
    // the ‖r‖ query). cd may be null ⇒ a pure penalty tie (λ_tie≡0, no global accumulation).
    void addMortarTieForce(const double D[4][4], const double M[4][4],
                           class LadrunoContactDomain *cd);
};

#endif
