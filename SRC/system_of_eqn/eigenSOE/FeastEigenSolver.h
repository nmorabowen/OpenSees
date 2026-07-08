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

// Ladruno: FEAST band-targeted eigensolver — solver twin (ADR 43, P1). See
// Ladruno_implementation/43_ladruno_feast_eigensolver_adr.md.
//
// Drives MKL's Extended Eigensolver (dfeast_scsrgv) on the CSR K,M owned by
// FeastEigenSOE to return every mode whose eigenvalue lambda (= omega^2) lies
// in the band [emin, emax]. The contour/quadrature machinery lives entirely
// inside MKL FEAST; this class only marshals the arrays, runs the subspace-
// saturation enlargement loop, and stores the accepted pairs ascending.
//
// P3c-serial (-rci): alternatively drives the FEAST RCI (dfeast_srci) with
// LadrunoBlockZKernel as the inner complex-shift solve — same eigenpairs as
// the driver path; the orchestration seam the MPI rung distributes.
//
// MKL provides the FEAST kernels; on this fork's build matrix MKL is present
// only on the Windows/oneAPI build, so the solve body is guarded on _WIN32.
// Every other platform compiles a stub that reports "FEAST requires MKL".

#ifndef FeastEigenSolver_h
#define FeastEigenSolver_h

#include <EigenSolver.h>
#include <FeastEigenSOE.h>

class FeastEigenSolver : public EigenSolver
{
  public:
    FeastEigenSolver();
    virtual ~FeastEigenSolver();

    virtual int solve(int numModes, bool generalized, bool findSmallest = true);
    virtual int setSize(void);
    virtual int setEigenSOE(FeastEigenSOE &theSOE);

    virtual const Vector &getEigenvector(int mode);
    virtual double getEigenvalue(int mode);

    // number of accepted (reported) modes from the last solve — the band
    // defines the count, so the command layer reconciles the domain's
    // eigenvalue list to this after the analysis-driven solve
    int getNumFoundModes(void) const {return numModes;}

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

  protected:

  private:
    FeastEigenSOE *theSOE;
    int numModes;         // accepted & reported modes (ascending)
    double *eigenvalues;  // length numModes
    double *eigenvectors; // n * numModes, column-major (mode k at [k*n])
    int valueCapacity;    // allocated length of eigenvalues
    int vectorCapacity;   // allocated length of eigenvectors
    Vector *eigenV;       // persistent per-mode return vector (size n)
};

#endif
