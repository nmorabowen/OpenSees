/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 07/2026
//
// LadrunoUP — the unified Biot u-p (solid displacement + pore pressure)
// saturated-porous continuum element of the Ladruno fork (ADR 71, P1).
// ONE element class covers the whole geometry family through the stateless
// shape providers of LadrunoUPShapes.h — (ndm, nNodes) selects the provider:
//   (2,3) T3 · (2,4) Q4 · (3,8) H8 · (2,6) Bézier T6 · (3,10) Bézier Tet10
// P1 gates the Q4 lane end-to-end; T3/H8 parse and run (gates at P2); the
// Bézier Taylor–Hood lanes are guarded off until P3 (heterogeneous-ndf
// plumbing). All coupled-block math (Q/H/S/H̃/f_seep) lives in the shared
// pure kernel LadrunoUPKernel.h — this class does material handling, block
// placement into the OpenSees element APIs, loads/responses/parameters.
//
// THE HONEST-p CONTRACT (ADR 71 §3.2 — the load-bearing divergence from the
// upstream UP-ucsd/UWelements families): the nodal pressure DOF *is* p (not
// ∫p·dt). Blocks split by the time-derivative they multiply:
//
//   getTangentStiff()            [ K  , −Q ; 0 , H ]
//   getInitialStiff()            [ K₀ , −Q ; 0 , H ]   (cached; p-rows NONZERO)
//   getDamp()                    [ C_Ray , 0 ; Qᵀ , S + H̃ ]
//   getMass()                    [ M , 0 ; 0 , 0 ]     (consistent | -lumped HRZ row-sum)
//   getResistingForce()          [ ∫BᵀσʹdV − Q·p − body ; H·p − f_seep ]  (NO rate terms)
//   getResistingForceIncInertia  = resisting + M·ü + getDamp()·v
//                                  (rates from trial nodal vel/accel; the p-slot
//                                   "velocity" is ṗ — integrators never assemble
//                                   damp·vel themselves, the element must)
//   addInertiaLoadToUnbalance    −M·(R·üg), p-slots zeroed
//
// Rayleigh is SOLID-ONLY by contract: this class overrides getDamp,
// getResistingForceIncInertia, setRayleighDampingFactors and shadows
// getRayleighDampingForces — C_Ray = αM·M_s + βK·K_s(current) + βK0·K_s(initial)
// + βKc·K_s(committed), u-u block only, with the element keeping its OWN
// committed/initial solid-K copies (never the base-class helpers, which would
// compose βK·getTangentStiff() and inject −Q/H into the damping). Porting
// reference: SSPquadUP.cpp:437–449 — noting it folds betaK0/betaKc onto the
// CURRENT solid K, which we fix here (true K₀/Kc copies).
//
// The effective-stress seam: update() maps trial u → strains → setTrialStrain
// ONLY. p never enters the constitutive update; it acts through Q in assembly.
//
// The price of honesty: the effective transient tangent is UNSYMMETRIC (−Q and
// Qᵀ live in different matrices) ⇒ general solver required (UmfPack / SuperLU /
// FullGeneral / BandGeneral; MUMPS SYM=0 under MPI). The default ProfileSPD
// silently drops one Q block (parser notices at creation; documented quirk).
//
// See Ladruno_implementation/71_ladruno_up_family_adr.md and
// Ladruno_implementation/_adr71_implementation_plan.md (WP1.A pins).

#ifndef LadrunoUP_h
#define LadrunoUP_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <vector>

class Node;
class NDMaterial;
class Response;

class LadrunoUP : public Element
{
  public:
    // shape provider dispatch (stateless P0 structs; no virtual provider class)
    enum ShapeKind { SHAPE_NONE = 0, SHAPE_T3, SHAPE_Q4, SHAPE_H8,
                     SHAPE_BT6, SHAPE_BTET10 };

    // WP1.A pinned ctor (plan §"WP1.A pins", FROZEN 2026-07-10)
    LadrunoUP(int tag, int ndm, const ID &nodeTags,        // 3|4|8 nodes at P1
              NDMaterial &theMaterial,                     // copied per GP
              double thickness,                            // 2D only; 1.0 in 3D
              double Kf, double poro, double rhoF,
              const Vector &perm,                          // ndm entries, k_hydraulic/gammaW
              double biotAlpha, double Ks,                 // Ks<=0 => infinite grains
              const Vector &body, const Vector &fluidBody, // ndm entries each
              int formulation,                             // 0=std, 1=bbar
              int pOrder,                                  // 0=equal (P1); 1=TH reserved P3
              bool lumped,
              int stabMode, double stabValue,              // 0=off,1=auto(alpha0),2=manual(alpha)
              bool dynSeepage);                            // default true
    LadrunoUP();                                           // broker
    ~LadrunoUP();

    const char *getClassType(void) const { return "LadrunoUP"; }

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);
    double getCharacteristicLength(void);   // crack-band lch = sqrt(A) | cbrt(V)

    // state
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // stiffness / damping / mass / residual — the six §3.2 contract APIs
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);
    const Matrix &getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    // solid-only Rayleigh bookkeeping (contract obligation 2, ADR §3.2)
    int setRayleighDampingFactors(double alphaM, double betaK,
                                  double betaK0, double betaKc);

    // parallel / IO — STUBS at P1 wave 1; bodies land in wave 2 (WP1.D)
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

    // updateMaterialStage transport + xPerm/yPerm/zPerm (contract obligation 3)
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

  protected:
    // shadows the (non-virtual) base helper: SOLID-ONLY Rayleigh forces,
    // never βK·getTangentStiff() (which carries −Q/H under this contract).
    const Vector &getRayleighDampingForces(void);

  private:
    // ---- construction / reconstruction (ctor + recvSelf share this) --------
    // provider dispatch + dofMap + all size-derived members + per-instance
    // algebra/cache allocation, driven purely from {ndm_, nNodes_, pOrder_} and
    // the physical params already set as members. Called by the real ctor AND
    // by recvSelf so the two paths can never drift (ADR 71 §3.5 anti-drift).
    // On an illegal (ndm_, nNodes_) it leaves shapeKind_ == SHAPE_NONE (caller
    // guards). Does NOT touch connectedExternalNodes_, the materials, or the
    // ctor-object params (perm/body/fluidBody) — those are set by the caller.
    void configureSizing(void);

    // ---- shape-provider dispatch (switch on shapeKind_) --------------------
    int  shapeNGP(void) const;
    void shapeEval(int gp, double *Nu, double *dNu, double *detJ) const;
    void shapeEvalP(int gp, bool taylorHood, double *Np, double *dNp) const;
    double shapeElemSizeH(void) const;
    const int *shapePCarrierTH(void) const;
    double shapeGPWeight(int gp) const;

    // ---- assembly helpers ---------------------------------------------------
    void formB(int gp, Matrix &B) const;      // solid B (mean-dilatation B-bar if formulation==1)
    void buildSolidK(Matrix &Ks, bool initial);         // ∑ BᵀDB dv (current|initial tangent)
    void buildSolidMass(Matrix &Ms) const;    // consistent | HRZ row-sum (lumped_)
    void buildCRaySolid(Matrix &C);           // αM·M_s + βK·K_s + βK0·K_s0 + βKc·K_sc
    void buildStaticBlocks(void);             // Q/H/S/H̃(unit-α) from kernel, geometry-fixed
    double stabAlphaValue(void);              // lazy auto-α (McGann) with cache + dirty flag
    void gpAccel(int gp, double *acc) const;  // trial solid acceleration interpolated at GP
    double gpPressure(int gp) const;          // Np · p_nodal (trial)

    // ---- topology / configuration ------------------------------------------
    int ndm_;                        // 2 | 3
    int nNodes_;                     // 3 | 4 | 8 (| 6 | 10 guarded to P3)
    int shapeKind_;                  // ShapeKind
    int nGP_;
    int nU_;                         // nNodes_ * ndm_  (local solid DOFs)
    int nP_;                         // # pressure carriers (= nNodes_ equal-order)
    int nStr_;                       // 3 (2D plane strain) | 6 (3D)
    int nDof_;                       // total element DOFs (kernel dofMap)

    ID connectedExternalNodes_;
    Node *theNodes_[10];
    NDMaterial **theMaterials_;      // nGP_ per-GP copies

    // ---- physical parameters (ctor) -----------------------------------------
    double thickness_;
    double Kf_, poro_, rhoF_;
    double kbar_[3];                 // diagonal Darcy k_hydraulic/γ_w
    double biotAlpha_, Ks_;          // grain bulk (Ks<=0 => infinite)
    double b_[3], fluidBody_[3];     // ctor body / fluid-body accelerations
    int formulation_;                // 0=std, 1=bbar
    int pOrder_;                     // 0=equal, 1=TH (P3)
    bool lumpedMass_;
    int stabMode_;                   // 0=off, 1=auto, 2=manual
    double stabValue_;               // alpha0 (auto) | alpha (manual)
    bool dynSeepage_;
    double oneOverQbar_;             // storage coefficient (kernel)

    // ---- loads --------------------------------------------------------------
    int applyLoad_;
    double appliedB_[3], appliedFluidBody_[3];
    Vector Qload_;                   // applied nodal loads (inertia unbalance)

    // ---- per-instance matrices/vectors (NO class statics — ADR-40) ----------
    Matrix K_;                       // tangent/initial scratch (nDof×nDof)
    Matrix damp_;                    // damping (nDof×nDof)
    Matrix mass_;                    // mass (nDof×nDof)
    Vector resid_;                   // resisting force (nDof)
    Vector rayForce_;                // Rayleigh force scratch (nDof)

    // solid-K bookkeeping for solid-only Rayleigh (ADR §3.2 obligation 2)
    Matrix KinitSolid_;              // initial solid K (βK0); lazily built
    bool   kInitSolidValid_;
    Matrix KcSolid_;                 // committed solid K (βKc); refreshed in commitState
    bool   kcSolidValid_;

    // caches with dirty flags (updateMaterialStage / perm setParameter)
    Matrix K0_;                      // full [K₀,−Q;0,H]
    bool   k0Valid_;
    bool   k0Dirty_;
    double stabAlpha_;
    bool   stabDirty_;

    // ---- geometry caches (linear geometry — fixed after setDomain) ----------
    bool geomValid_;
    std::vector<double> xy_;         // nodal coords, xy_[a*ndm+i]
    std::vector<double> NuAll_;      // nGP*nNodes
    std::vector<double> dNuAll_;     // nGP*nNodes*ndm (cartesian)
    std::vector<double> NpAll_;      // nGP*nP
    std::vector<double> dNpAll_;     // nGP*nP*ndm
    std::vector<double> dv_;         // nGP: w·detJ(·thick in 2D)
    std::vector<double> dNuBar_;     // nNodes*ndm element-mean gradients (B-bar)
    double elemH_;                   // provider elemSizeH (largest edge)
    double charLen_;                 // sqrt(A) | cbrt(V)
    std::vector<int> uOff_, pOff_;   // kernel dofMap slots
    std::vector<int> carrierNodes_;  // node index of each pressure carrier (len nP_)

    // geometry-fixed kernel blocks (row-major, kernel layout)
    std::vector<double> Qblk_;       // nU × nP   (B-bar-consistent under bbar)
    std::vector<double> Hblk_;       // nP × nP   (rebuilt when perm changes)
    std::vector<double> Sblk_;       // nP × nP
    std::vector<double> HtUnit_;     // nP × nP   unit-α stabilization Laplacian

    // small per-instance assembly scratch (sized in setDomain)
    Matrix Bscratch_;                // nStr × nU
    Matrix KsolScratch_;             // nU × nU
    Matrix CRayScratch_;             // nU × nU
    Matrix MsolScratch_;             // nU × nU
    Vector epsScratch_;              // nStr
    Vector uSolScratch_;             // nU
    Vector fuScratch_;               // nU
    Vector pScratch_;                // nP
};

#endif
