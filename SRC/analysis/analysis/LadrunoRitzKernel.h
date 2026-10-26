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

// Ladruno WP-213: shared block Krylov-Ritz kernel on (K, M) (ADR 213 D1).
// See Ladruno_implementation/213_ladruno_modal_basis_completeness_adr.md.
//
//   One kernel serves `residualModes` (WP-213, missing-mass static-correction
//   pseudo-modes) and `ritzVectors` (PR 3, load-dependent Ritz vectors):
//
//     factor K once                     (the analysis' LinearSOE; the solves
//                                        after the first reuse the factor)
//     X0 = K^-1 F                       (q starting loads; or X0 = F for shape
//                                        starts, Options::solveStart = false)
//     for each block b:
//        deflate Xb against Phi and X<b (twice, M-inner product)
//        M-orthonormalize Xb             (modified Gram-Schmidt, twice; drop a
//                                        column whose M-norm fell below
//                                        dropTol * its norm right after the solve)
//        X(b+1) = K^-1 (M Xb)
//     K~ = X^T K X, eig(K~) -> (lambda, Z) (dense LAPACK dsyevx)
//     Psi = X Z                          (M-orthonormal, K-orthogonal)
//
//   Wilson, Yuan & Dickens (1982) in block form; Joo, Wilson & Leger (1989)
//   for the drop / breakdown rule.
//
//   K is what `eigen` assembles: every FE_Element's addKtToTang(1.0) at the
//   current state (so a gravity + loadConst + PDelta state gives K - Kg), plus
//   the constraint FEs of Penalty / Lagrange. M products go through
//   FE_Element::getM_Force / DOF_Group::getM_Force in equation numbering, the
//   inertia-load path, so Transformation (T^T M T), Penalty and Lagrange
//   (multiplier rows carry no mass) are all handled by the handler's own
//   objects. No classTag: the kernel is not a domain object.
//
//   Equation-numbered vectors include Lagrange multiplier rows (zero mass).
//   Every M-inner product ignores them; X^T K X does not see them either
//   because the constraint rows of a solve are satisfied (C u = 0).

#ifndef LadrunoRitzKernel_h
#define LadrunoRitzKernel_h

#include <vector>
#include <string>
#include <Vector.h>

class AnalysisModel;
class LinearSOE;
class Domain;
class Node;

// ---------------------------------------------------------------------------
// Mode-kind registry (ADR 213 D2): what each mode of the domain's CURRENT mode
// set is. Kept on the Ladruno side, keyed by Domain*, so no vanilla Domain
// member is needed. A record is only valid while the domain's eigenvalue
// vector is bit-identical to the snapshot taken when the record was written:
// a new `eigen` (or `wipe`, FEAST trim, ...) replaces that vector, so the
// record goes stale and every mode reads as EIGEN again. That is how "a new
// eigen clears the tags" without touching the eigen path.
// ---------------------------------------------------------------------------
namespace LadrunoModeKinds
{
    enum Kind { EIGEN = 0, RESIDUAL = 1, RITZ = 2 };

    // kinds of the domain's current modes (size = number of eigenvalues);
    // all EIGEN when there is no valid record
    std::vector<int> get(const Domain* theDomain);
    // record kinds for the domain's current mode set (snapshots its
    // eigenvalues); saPeriods (optional, same length) holds a spectral-ordinate
    // period override per mode, 0 = none
    void set(const Domain* theDomain, const std::vector<int>& kinds,
             const std::vector<double>& saPeriods = std::vector<double>(),
             const std::vector<int>& dirMasks = std::vector<int>());
    // Sa-period override of mode mode0 (0-based) of the domain's current set,
    // or 0 when there is none or the record is stale. responseSpectrumAnalysis
    // reads Sa at this period while lambda (u = phi Gamma Sa / lambda) and the
    // CQC omegas stay the mode's own (ADR 213 D3b)
    double saPeriod(const Domain* theDomain, int mode0);
    std::vector<double> saPeriods(const Domain* theDomain);
    // per mode: bit (d-1) set when residualModes -dir d produced it (0 for
    // eigen modes); lets -replace rebuild every direction it drops
    std::vector<int> dirMasks(const Domain* theDomain);
    void clear(const Domain* theDomain);
    // single-letter label for printing: E / R / Z (Ritz)
    char label(int kind);
    // Ladruno WP-215: true when some mode of the domain's current set is of
    // this kind (valid record only)
    bool has(const Domain* theDomain, int kind);
    // Ladruno WP-215: one-line WARNING on opserr when the current set holds
    // Ritz vectors (frequency-domain commands: a Ritz basis is tuned to a
    // load, not to a frequency band); returns has(theDomain, RITZ)
    bool warnRitzBasis(const Domain* theDomain, const char* cmd);

    // Ladruno WP-214 (ADR 213 D6): decoration of the modalProperties report.
    // Active only when the domain's record is valid, its eigenvalue snapshot
    // equals `lambdas` bit for bit and some mode is not EIGEN; otherwise every
    // string is empty (and kinds too), so the report stays the stock one.
    struct Report
    {
        std::vector<int> kinds;         // per mode; empty when inactive
        std::vector<std::string> tags;  // per mode: "" or "  # R" / "  # Z"
        std::string note;               // "# NOTE: ..." lines; "" when inactive
        bool active() const { return !kinds.empty(); }
        // tag of mode mode0 (0-based), "" when inactive / out of range
        const char* tag(int mode0) const;
    };
    Report report(const Domain* theDomain, const Vector& lambdas);
}

class LadrunoRitzKernel
{
public:
    struct Options
    {
        double dropTol = 1.0e-10;   // relative M-norm drop tolerance
        int    maxBlocks = 1;       // blocks per starting vector (>= 1)
        int    maxVectors = 0;      // stop once this many columns survive (0 = no cap)
        // true: F are loads and block 0 is K^-1 F (load-dependent Ritz, LDR).
        // false: F are shapes and block 0 is F itself, block 1 = K^-1 M F, ...
        // (residualModes starts from the rigid shape iota_d, so the deflated
        // residual shape AND its static response are both in the span)
        bool   solveStart = true;
        // optional per-start cap on blocks (CSI `ncyc`); empty = maxBlocks for all
        std::vector<int> cyclesPerStart;
    };

    struct Result
    {
        std::vector<Vector> psi;      // Ritz vectors, equation numbering, M-orthonormal
        std::vector<double> lambda;   // Ritz values, ascending, all > 0
        // per starting vector: how many of its Krylov columns survived / were dropped
        std::vector<int> keptPerStart;
        std::vector<int> droppedPerStart;
        int numDropped = 0;
        // max normwise backward error ||K x - b|| / (|| |K||x| || + ||b||) over
        // the first batch of solves (WP-215, #1025: ||K x - b|| / ||b|| is not
        // scale-free under Penalty, where || |K||x| || ~ alpha ||x||)
        double maxSolveResidual = 0.0;
        // per starting vector: M-norm of its block-0 column after deflation
        // over its M-norm before (0 when dropped); squared, the share of that
        // start's mass the existing basis misses (shape starts)
        std::vector<double> firstBlockRelNorm;
    };

    LadrunoRitzKernel(AnalysisModel& theModel, LinearSOE& theSOE);

    int numEqn() const { return m_neq; }

    // ---- products in equation numbering ----------------------------------
    // y = M x (FE_Element / DOF_Group getM_Force, constraint-aware)
    void applyM(const Vector& x, Vector& y);
    // KX[j] = K X[j] for all columns in one sweep over the FE_Elements.
    // absQuad (optional): absQuad[j] = sum_e |x_e|^T |K_e| |x_e|, the scale of
    // the round-off in x^T K x (huge next to x^T K x under Penalty).
    // absKX (optional, WP-215): absKX[j] = sum_e |K_e| |x_e| in equation
    // numbering, the scale of the round-off in K x (normwise backward error)
    void applyK(const std::vector<Vector>& X, std::vector<Vector>& KX,
                std::vector<double>* absQuad = 0,
                std::vector<Vector>* absKX = 0);
    // M-inner product with a precomputed M y
    static double dot(const Vector& a, const Vector& b);
    // Ladruno WP-215: diag(M) in equation numbering (element and nodal mass
    // matrices, the eigen assembly); offDiagonal is set when any assembled
    // mass matrix has a nonzero off-diagonal entry (consistent mass, or T^T M T
    // coupling under Transformation)
    void massDiagonal(Vector& d, bool& offDiagonal);

    // ---- domain <-> equation numbering ----------------------------------
    // rigid translation along global direction dir (1-based, <= ndm), zero at
    // every DOF a domain SP constraint fixes (those move with the ground; under
    // Penalty / Lagrange they still carry an equation)
    int rigidInfluence(int dir, Vector& iota);
    // Ladruno WP-215: the reference nodal loads of a load pattern assembled to
    // equation numbering, f = sum_groups T^T p_node (TransformationDOF_Group's
    // own T, so a load on a constrained DOF moves to the retained node; loads
    // on SP-fixed DOFs are dropped). numElementalLoads receives the number of
    // elemental loads the pattern holds (not assembled).
    int patternLoad(int patternTag, Vector& f, int& numElementalLoads);
    // the domain's current modes (node eigenvectors) as equation vectors.
    // Under Transformation the retained DOFs carry them exactly (constrained
    // nodes are recovered from T on the way back in).
    int domainModes(std::vector<Vector>& phi);
    // install a whole mode set into the domain exactly as eigen does:
    // setNumEigenvectors / setEigenvector / setEigenvalues
    int install(const std::vector<Vector>& modes, const std::vector<double>& lambdas);

    // ---- the kernel -------------------------------------------------------
    // Assemble K into the analysis' LinearSOE (the factorization happens on the
    // first solve and is reused by the later ones).
    int formTangent();
    // x = K^-1 b on the factor formed by formTangent()
    int solve(const Vector& b, Vector& x);

    // Block Krylov-Ritz on the starting vectors F (load vectors, equation
    // numbering), deflated against the optional basis Phi (need not be
    // normalized; it is M-orthonormalized internally). Forms K itself.
    int run(const std::vector<Vector>& F,
            const std::vector<Vector>* Phi,
            const Options& opt,
            Result& out);

    const std::string& error() const { return m_err; }

private:
    // map a per-node DOF vector (ndf values) into equation numbering
    template <class Getter> void nodeToEq(Getter g, Vector& out);

    AnalysisModel* m_model;
    LinearSOE*     m_soe;
    Domain*        m_domain;
    int            m_neq;
    std::string    m_err;
};

#endif
