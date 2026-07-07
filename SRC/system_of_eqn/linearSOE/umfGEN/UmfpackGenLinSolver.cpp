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
                                                                        
// $Revision: 1.4 $
// $Date: 2009-05-20 17:30:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/98
//
// Description: This file contains the class definition for
// UmfpackGenLinSolver. It solves the UmfpackGenLinSOEobject by calling
// UMFPACK5.7.1 routines.
//
// What: "@(#) UmfpackGenLinSolver.C, revA"

#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>
#include <math.h>
#include <cstring>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <profiler/ProfilerMacros.h>  // Ladruno (ADR-40 rank 8/10): factor-vs-trisolve split

void* OPS_UmfpackGenLinSolver()
{
    bool useLongIndices = false;
    // Ladruno (ADR-40 rank 2): expose UMFPACK strategy + pivot tolerance.
    // Default strategy AUTO (was hardcoded SYMMETRIC — mis-orders unsymmetric
    // tangents, e.g. LadrunoConcrete3D); default pivotTol keeps the legacy 1.0
    // (maximal-threshold partial pivoting). Legacy behavior is exactly
    //   system UmfPack -strategy symmetric -pivotTol 1.0
    int strategy = UMFPACK_STRATEGY_AUTO;
    double pivotTol = 1.0;
    // Ladruno (ADR-40 rank 8/10): persist the Numeric LU across solves by
    // default; -noNumericPersist restores the legacy free-every-solve behavior.
    bool persist = true;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (opt == nullptr) {
            opt = "";
        }
        if (strcmp(opt, "useLongIndices") == 0 ||
            strcmp(opt, "-useLongIndices") == 0) {
            useLongIndices = true;
        } else if (strcmp(opt, "-noNumericPersist") == 0) {
            persist = false;   // Ladruno (ADR-40 rank 8/10)
        } else if (strcmp(opt, "-strategy") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING UmfpackGenLinSolver: -strategy requires auto|symmetric|unsymmetric\n";
                return nullptr;
            }
            const char *s = OPS_GetString();
            if (s == nullptr) s = "";
            if (strcmp(s, "auto") == 0) {
                strategy = UMFPACK_STRATEGY_AUTO;
            } else if (strcmp(s, "symmetric") == 0) {
                strategy = UMFPACK_STRATEGY_SYMMETRIC;
            } else if (strcmp(s, "unsymmetric") == 0) {
                strategy = UMFPACK_STRATEGY_UNSYMMETRIC;
            } else {
                opserr << "WARNING UmfpackGenLinSolver: unknown -strategy \"" << s
                       << "\" (expected auto|symmetric|unsymmetric)\n";
                return nullptr;
            }
        } else if (strcmp(opt, "-pivotTol") == 0) {
            int numData = 1;
            if (OPS_GetDoubleInput(&numData, &pivotTol) != 0) {
                opserr << "WARNING UmfpackGenLinSolver: -pivotTol requires a float\n";
                return nullptr;
            }
        } else {
            opserr << "WARNING UmfpackGenLinSolver: unknown option \"" << opt
                   << "\" (expected useLongIndices, -noNumericPersist, -strategy, or -pivotTol)\n";
            return nullptr;
        }
    }

    UmfpackGenLinSolver *theSolver = new UmfpackGenLinSolver(useLongIndices, strategy, pivotTol, persist);
    return new UmfpackGenLinSOE(*theSolver);
}

UmfpackGenLinSolver::
UmfpackGenLinSolver(bool useLongIndices, int strategy, double pivotTol, bool persist)
    : LinearSOESolver(SOLVER_TAGS_UmfpackGenLinSolver),
      useLongIndices(useLongIndices),
      strategy(strategy),      // Ladruno (ADR-40 rank 2)
      pivotTol(pivotTol),      // Ladruno (ADR-40 rank 2)
      persist(persist),        // Ladruno (ADR-40 rank 8/10)
      Symbolic(0),
      Numeric(0),              // Ladruno (ADR-40 rank 8/10)
      theSOE(0)
{
}

UmfpackGenLinSolver::~UmfpackGenLinSolver()
{
    this->freeNumeric();   // Ladruno (ADR-40 rank 8/10)
    if (Symbolic != 0) {
        if (useLongIndices) {
#ifdef _UMFPACK_DLONG
            umfpack_dl_free_symbolic(&Symbolic);
#endif
        } else {
            umfpack_di_free_symbolic(&Symbolic);
        }
    }
}

// Ladruno (ADR-40 rank 8/10): free the persisted Numeric LU (both index modes).
void
UmfpackGenLinSolver::freeNumeric(void)
{
    if (Numeric != 0) {
        if (useLongIndices) {
#ifdef _UMFPACK_DLONG
            umfpack_dl_free_numeric(&Numeric);
#endif
        } else {
            umfpack_di_free_numeric(&Numeric);
        }
        Numeric = 0;
    }
}

void
UmfpackGenLinSolver::syncIndexBuffers(void)
{
    if (!useLongIndices) return;

    const int n = theSOE->X.Size();
    const std::size_t nnz = theSOE->Ai.size();
    Ap64.resize(static_cast<std::size_t>(n) + 1u);
    Ai64.resize(nnz);
    for (std::size_t i = 0; i < Ap64.size(); ++i) {
        Ap64[i] = static_cast<SuiteSparse_long>(theSOE->Ap[i]);
    }
    for (std::size_t k = 0; k < nnz; ++k) {
        Ai64[k] = static_cast<SuiteSparse_long>(theSOE->Ai[k]);
    }
}

int
UmfpackGenLinSolver::solve(void)
{
    int n = theSOE->X.Size();
    int nnz = (int)theSOE->Ai.size();
    if (n == 0 || nnz==0) return 0;

    double* Ax = theSOE->Ax.data();
    double* X = &(theSOE->X(0));
    double* B = &(theSOE->B(0));

    // check if symbolic is done
    if (Symbolic == 0) {
        opserr<<"WARNING: setSize has not been called -- Umfpackgenlinsolver::solve\n";
        return -1;
    }

    // Ladruno (ADR-40 rank 8/10): reuse the persisted Numeric LU while the
    // assembled matrix is unchanged. theSOE->factored is the standard direct-
    // solver invalidation flag (reset by zeroA/setSize, i.e. every tangent
    // reassembly). We refactor iff persistence is off, the matrix was
    // reassembled, or we hold no valid Numeric (fresh/swapped-in solver).
    // The Numeric!=0 guard makes a spuriously-true factored flag harmless.
    const bool needFactor = (!persist) || (!theSOE->factored) || (Numeric == 0);

    if (useLongIndices) {
#ifdef _UMFPACK_DLONG
        SuiteSparse_long *Ap = Ap64.data();
        SuiteSparse_long *Ai = Ai64.data();
        SuiteSparse_long status;

        if (needFactor) {
            this->freeNumeric();   // drop any stale factorization first
            { OPS_PROFILE_SCOPE("soe.factor");   // Ladruno (ADR-40 rank 8/10)
            status =
                umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
            }
            if (status != UMFPACK_OK) {
                opserr << "WARNING: numeric analysis returns "
                       << static_cast<int>(status)
                       << " -- Umfpackgenlinsolver::solve\n";
                this->freeNumeric();
                theSOE->factored = false;
                return -1;
            }
            theSOE->factored = true;
        }

        // solve (always) using the current/persisted factorization
        { OPS_PROFILE_SCOPE("soe.trisolve");   // Ladruno (ADR-40 rank 8/10)
        status = umfpack_dl_solve(UMFPACK_A, Ap, Ai, Ax, X, B, Numeric, Control,
                                  Info);
        }

        if (status != UMFPACK_OK) {
            opserr << "WARNING: solving returns " << static_cast<int>(status)
                   << " -- Umfpackgenlinsolver::solve\n";
            this->freeNumeric();
            theSOE->factored = false;
            return -1;
        }

        // legacy free-every-solve behavior when persistence is disabled
        if (!persist) { this->freeNumeric(); theSOE->factored = false; }
#endif
    } else {
        int *Ap = theSOE->Ap.data();
        int *Ai = theSOE->Ai.data();
        int status;

        if (needFactor) {
            this->freeNumeric();   // drop any stale factorization first
            { OPS_PROFILE_SCOPE("soe.factor");   // Ladruno (ADR-40 rank 8/10)
            status =
                umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
            }
            if (status != UMFPACK_OK) {
                opserr << "WARNING: numeric analysis returns " << status
                       << " -- Umfpackgenlinsolver::solve\n";
                this->freeNumeric();
                theSOE->factored = false;
                return -1;
            }
            theSOE->factored = true;
        }

        // solve (always) using the current/persisted factorization
        { OPS_PROFILE_SCOPE("soe.trisolve");   // Ladruno (ADR-40 rank 8/10)
        status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, X, B, Numeric, Control,
                                  Info);
        }

        if (status != UMFPACK_OK) {
            opserr << "WARNING: solving returns " << status
                   << " -- Umfpackgenlinsolver::solve\n";
            this->freeNumeric();
            theSOE->factored = false;
            return -1;
        }

        // legacy free-every-solve behavior when persistence is disabled
        if (!persist) { this->freeNumeric(); theSOE->factored = false; }
    }

    return 0;
}

int
UmfpackGenLinSolver::setSize()
{
    // Ladruno (ADR-40 rank 8/10): the sparsity structure (and thus Symbolic) is
    // about to be rebuilt -> any persisted Numeric is stale. Free it here; the
    // SOE has already reset factored=false.
    this->freeNumeric();

    int n = theSOE->X.Size();
    int nnz = (int)theSOE->Ai.size();
    if (n == 0 || nnz == 0) {
        Ap64.clear();
        Ai64.clear();
        return 0;
    }

    if (useLongIndices) {
#ifdef _UMFPACK_DLONG
        // set default control parameters
        umfpack_dl_defaults(Control);
        Control[UMFPACK_PIVOT_TOLERANCE] = pivotTol;   // Ladruno (ADR-40 rank 2): was hardcoded 1.0
        Control[UMFPACK_STRATEGY] = strategy;          // Ladruno (ADR-40 rank 2): was hardcoded SYMMETRIC

        syncIndexBuffers();
        SuiteSparse_long *Ap = Ap64.data();
        SuiteSparse_long *Ai = Ai64.data();
        double *Ax = theSOE->Ax.data();

        // symbolic analysis
        if (Symbolic != 0) {
            umfpack_dl_free_symbolic(&Symbolic);
        }
        OPS_PROFILE_SCOPE("soe.symbolic");   // Ladruno (ADR-40 rank 8/10)
        SuiteSparse_long status = umfpack_dl_symbolic(
            (SuiteSparse_long)n, (SuiteSparse_long)n, Ap, Ai, Ax, &Symbolic,
            Control, Info);

        // check error
        if (status != UMFPACK_OK) {
            opserr << "WARNING: symbolic analysis returns "
                   << static_cast<int>(status)
                   << " -- Umfpackgenlinsolver::setsize\n";
            Symbolic = 0;
            return -1;
        }
#endif
    } else {
        // set default control parameters
        umfpack_di_defaults(Control);
        Control[UMFPACK_PIVOT_TOLERANCE] = pivotTol;   // Ladruno (ADR-40 rank 2): was hardcoded 1.0
        Control[UMFPACK_STRATEGY] = strategy;          // Ladruno (ADR-40 rank 2): was hardcoded SYMMETRIC

        int *Ap = theSOE->Ap.data();
        int *Ai = theSOE->Ai.data();
        double *Ax = theSOE->Ax.data();

        // symbolic analysis
        if (Symbolic != 0) {
            umfpack_di_free_symbolic(&Symbolic);
        }
        OPS_PROFILE_SCOPE("soe.symbolic");   // Ladruno (ADR-40 rank 8/10)
        int status =
            umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, Control, Info);
        
        
        // check error
        if (status != UMFPACK_OK) {
            opserr << "WARNING: symbolic analysis returns " << status
                   << " -- Umfpackgenlinsolver::setsize\n";
            Symbolic = 0;
            return -1;
        }
    }

    return 0;
}

int
UmfpackGenLinSolver::setLinearSOE(UmfpackGenLinSOE &theLinearSOE)
{
    theSOE = &theLinearSOE;
    return 0;
}

int
UmfpackGenLinSolver::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
UmfpackGenLinSolver::recvSelf(int ctag, 
                              Channel &theChannel,
                              FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}
