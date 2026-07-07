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

// Ladruno: FEAST band-targeted eigensolver SOE (ADR 43, P1). See
// Ladruno_implementation/43_ladruno_feast_eigensolver_adr.md.
//
// A self-assembling EigenSOE for the real symmetric generalized pencil
//     K phi = lambda M phi,   K symmetric, M SPD,
// storing K and M in ONE shared CSR sparsity pattern (3-array, 1-based
// indices, FULL pattern -> MKL uplo = 'F') so the twin FeastEigenSolver can
// hand the arrays straight to MKL's Extended Eigensolver (dfeast_scsrgv),
// which returns every mode whose lambda falls in the requested band
// [emin, emax] (lambda units = (2*pi*f)^2 — the orchestrator converts Hz).
//
// The SOE owns its storage (like SymBandEigenSOE); the solver only reads it.
// The band / contour parameters live here so the parser can stash them on the
// SOE before solve(); the solver pulls them back at solve time.

#ifndef FeastEigenSOE_h
#define FeastEigenSOE_h

#include <EigenSOE.h>

class FeastEigenSolver;

class FeastEigenSOE : public EigenSOE
{
  public:
    // no AnalysisModel needed (ArpackSOE precedent) — lets the command layer
    // construct + configure the SOE BEFORE any analysis object exists
    FeastEigenSOE(FeastEigenSolver &theSolver);

    virtual ~FeastEigenSOE();

    virtual int getNumEqn(void) const;
    virtual int setSize(Graph &theGraph);

    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual int addM(const Matrix &, const ID &, double fact = 1.0);

    virtual void zeroA(void);
    virtual void zeroM(void);

    // --- band / contour parameters (lambda units = (2*pi*f)^2) ------------
    // The parser converts the user's Hz band before calling setBand().
    void setBand(double emin, double emax);
    void setSubspaceSize(int m0);   // 0 => auto (solver seeds a stochastic estimate)
    void setNumQuad(int nq);        // MKL fpm[1] = number of contour points
    void setMaxRefine(int maxRefine);
    void setTol(double eps);        // MKL fpm[2] = stopping exponent 10^-eps
    void setVerbose(bool verbose);  // MKL fpm[0] runtime print

    // found-mode count from the last solve (the band defines the count;
    // the command layer reconciles the domain's list to this)
    int getNumFoundModes(void) const;

    double getEmin(void) const;
    double getEmax(void) const;
    int    getSubspaceSize(void) const;
    int    getNumQuad(void) const;
    int    getMaxRefine(void) const;
    int    getTol(void) const;
    bool   getVerbose(void) const;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    friend class FeastEigenSolver;

  protected:

  private:
    // locate the value-array slot for (row,col) [0-based] within the shared
    // CSR pattern via binary search over the row's ascending column list;
    // returns the 0-based index into Kvals/Mvals, or -1 if not in pattern.
    int findSlot(int row, int col) const;

    int size;            // number of equations n
    int nnz;             // number of stored entries (active)
    int csrCapacity;     // allocated length of colInd/Kvals/Mvals
    int *rowPtr;         // n+1, 1-based row starts (MKL isa/isb)
    int *colInd;         // nnz, 1-based column indices, ascending per row
    double *Kvals;       // nnz, stiffness values (MKL sa)
    double *Mvals;       // nnz, mass values      (MKL sb)

    // FEAST / contour parameters
    double emin;
    double emax;
    int m0;              // subspace size (0 => auto)
    int nq;              // quadrature / contour points
    int maxRefine;       // max FEAST refinement loops
    int tolExp;          // stopping exponent (fpm[2] = 10^-tolExp)
    bool verbose;
};

#endif
