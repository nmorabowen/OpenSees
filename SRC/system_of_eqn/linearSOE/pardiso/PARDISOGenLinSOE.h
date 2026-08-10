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
                                                                        
                                                                       
#ifndef PARDISOGenLinSOE_h
#define PARDISOGenLinSOE_h

// Written: M. Salehi opensees.net@gmail.com
// website : http://opensees.net
// Created: 02/19
// Revision: A

#include <LinearSOE.h>
#include <Vector.h>

class PARDISOGenLinSolver;

// Ladruno ADR-75 P1d (2026-07): symmetric half-storage.
//
// `matType` selects what this SOE actually stores, and (through the solver,
// which reads it) the PARDISO `mtype`:
//
//   0  unsymmetric        full CSR                     -> mtype  11
//   1  symmetric pos.def. UPPER-triangle CSR           -> mtype   2
//   2  symmetric general  UPPER-triangle CSR           -> mtype  -2
//
// The numbering deliberately matches `system Mumps -matrixType 0|1|2` so the
// two solvers take the same verb. MKL PARDISO's symmetric modes require the
// upper triangle only, in row-major CSR, with ascending column indices and the
// diagonal ALWAYS present (even when zero) — all three invariants are
// established in setSize() below.
//
// WARNING (documented, deliberate): with matType != 0 only the col >= row half
// of each element matrix is read. If the assembled tangent is genuinely
// unsymmetric (contact, non-associated flow, follower loads, corotational
// transforms, LadrunoUP) this solves with the upper triangle reflected — it
// does not average. That is why the default stays 0.
//
// It DOES now detect it, though (ADR-75 P1d, adversarial review): addA already
// visits the col < row entries in order to discard them, so comparing each
// against its mirror costs one subtraction on data already in cache and turns
// the worst failure mode this feature has — a converged, plausible, WRONG
// answer — into a loud warning. See addA().

class PARDISOGenLinSOE : public LinearSOE
{
  public:
    PARDISOGenLinSOE(PARDISOGenLinSolver &theSolver);
    PARDISOGenLinSOE(PARDISOGenLinSolver &theSolver, int matType);

    ~PARDISOGenLinSOE();

    int getNumEqn(void) const;
    int getMatType(void) const;   // Ladruno ADR-75 P1d
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        
    
    void zeroA(void);
    void zeroB(void);
    
    const Vector &getX(void);
    const Vector &getB(void);    
    double normRHS(void);

    void setX(int loc, double value);        
    void setX(const Vector &x);        
    int setPARDISOGenLinSolver(PARDISOGenLinSolver &newSolver);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
	friend class PARDISOGenLinSolver;
  protected:
    
  private:
    int size;            // order of A
    int nnz;             // number of non-zeros in A
    double *A, *B, *X;   // 1d arrays containing coefficients of A, B and X
    int *colA, *rowStartA; // int arrays containing info about coefficients in A
    Vector *vectX;
    Vector *vectB;    
    int Asize, Bsize;    // size of the 1d array holding A
    bool factored;
    int matType;         // Ladruno ADR-75 P1d: 0 unsym / 1 SPD / 2 sym-general
    int asymWarned;      // Ladruno ADR-75 P1d: half-store asymmetry reported once
    int asymBudget;      // ...and how many more element matrices to check
    int asymPass;        // ...and how many tangent assemblies this pattern has
                         // seen, so the check can RE-ARM periodically and catch
                         // late-onset asymmetry (Ladruno TIMs item 2; see the
                         // ASYM_* constants and zeroA in the .cpp)
    int missWarned;      // Ladruno ADR-75 P1f: an addA entry with no CSR slot
                         // (stiffness silently discarded) reported once
};


#endif

