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

// Ladruno: complex / state-space modal analysis driver (ADR 46). See
// Ladruno_implementation/46_ladruno_complex_modal_adr.md.
//
//   P0 scope (this file): the REDUCED-PENCIL QZ KERNEL only — solve the
//   quadratic eigenproblem
//
//       (lambda^2 Mt + lambda Ct + Kt) z = 0,     Mt,Ct,Kt  p x p (small),
//
//   via the symmetric-block state-space linearization solved dense with
//   LAPACK dggev, and recover per-mode engineering quantities
//   (omega0, omegaD, zeta, complex reduced shape z). No Domain coupling —
//   the projection of assembled M, C, K onto the real modal basis Phi is
//   ADR 46 P1/P2; this kernel consumes the already-reduced matrices.
//
//   LINEARIZATION (ADR 46 §4.4, sign-corrected at P0):  with
//   w = {z ; lambda z}, the quadratic problem is equivalent to
//
//       B w = -lambda A w,   A = [Ct Mt; Mt 0],  B = [Kt 0; 0 -Mt].
//
//   We feed LAPACK the pair (-B, A):   (-B) w = lambda A w,  so dggev's
//   eigenvalue IS lambda directly. A is nonsingular whenever Mt is
//   (block anti-diagonal), so with a mass-normalized basis (Mt = I)
//   infinite eigenvalues (beta ~ 0) occur only for pathological input and
//   are detected + tagged, never divided through.
//
//   classTag 33019 (LadrunoComplexEigen) is RESERVED in
//   LEDGER_implementations.md; it enters SRC/classTags.h when ADR 46 P1
//   (the domain-coupled command) merges and the ledger row flips
//   RESERVED -> active.

#ifndef LadrunoComplexEigen_h
#define LadrunoComplexEigen_h

#include <vector>

class Matrix;

class LadrunoComplexEigen
{
  public:
    enum ModeKind {
        UNDERDAMPED = 0,  // conjugate pair, reported once (+Im representative)
        OVERDAMPED  = 1,  // real non-oscillatory root, reported per root
        RIGID       = 2   // lambda ~ 0, or infinite eigenvalue (beta ~ 0)
    };

    struct ComplexMode {
        double omega0;    // undamped natural circular frequency  |lambda|
        double omegaD;    // damped circular frequency  Im(lambda) (>= 0)
        double zeta;      // modal damping ratio  -Re(lambda)/|lambda|
        double lambdaRe;  // Re(lambda)
        double lambdaIm;  // Im(lambda) of the reported (+Im) representative
        ModeKind kind;
        double resid;     // ||(lambda^2 Mt + lambda Ct + Kt) z||_2, quality metric
        std::vector<double> zRe;  // reduced eigenvector z (top block), real part
        std::vector<double> zIm;  //                                    imag part
    };

    // Solve the reduced quadratic pencil. Mt, Ct, Kt must be p x p with the
    // same p. Underdamped conjugate pairs are collapsed to one physical mode
    // (the +Im root); overdamped/rigid real roots are each reported. Modes
    // are returned sorted by ascending omega0 (rigid first). z is
    // phase-normalized: the largest-|z_j| entry is made real-positive.
    // tol scales the beta ~ 0 (infinite-eigenvalue) and Im ~ 0
    // (non-oscillatory) classification tests.
    // Returns 0 on success; < 0 on dimension mismatch or LAPACK failure
    // (dggev INFO is reported via opserr).
    static int solveReducedPencil(const Matrix &Mt, const Matrix &Ct,
                                  const Matrix &Kt,
                                  std::vector<ComplexMode> &modes,
                                  double tol = 1.0e-8);
};

#endif
