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
** ****************************************************************** */

// Developed by: N. Mora Bowen, P. Palacios, J.A. Abell
//
// References:
//   Irons, B.M., "Applications of a theorem on eigenvalues to finite element
//     problems", University of Wales, Swansea, Report CR/132/70, 1970.
//   Hughes, T.J.R., Pister, K.S., Taylor, R.L., "Implicit-explicit finite elements
//     in nonlinear transient analysis", Computer Methods in Applied Mechanics and
//     Engineering, 17/18, 1979.
//   Flanagan, D.P., Belytschko, T., "Eigenvalues and stable time steps for the
//     uniform strain hexahedron and quadrilateral", Journal of Applied Mechanics,
//     51, 1984.
//
// Description: critical time step estimate from the per-element generalized
// eigenproblem K v = lambda M v. The maximum frequency of the assembled system
// does not exceed the maximum over the elements, so min_e 2/omega_e is a lower
// bound on the stability limit 2/omega_max of the central difference scheme. The
// bound is exact when an element's highest mode is also a mode of the assembled
// system, for example a free-free chain of equal lumped-mass bars.

#ifndef CriticalTimeStep_h
#define CriticalTimeStep_h

class AnalysisModel;
class Element;
class Matrix;
class Vector;

// Lumping of the element mass to the diagonal M of the per-element pencil.
enum class CTSLumping {
    RowSum,   // sum of each row on the diagonal; row sums of rotational DOFs of
              // beams and shells can vanish
    Diagonal  // diagonal of the consistent mass, the mass `system Diagonal` solves
              // with; positive, not mass-conserving
};

struct CTSResult {
    double damped_dt;      // limit including element Rayleigh damping; +inf if none
    double undamped_dt;    // 2/omega_max; +inf if none
    int    damped_tag;     // governing element (damped); -1 if none
    int    undamped_tag;   // governing element (undamped); -1 if none
    int    n_contributing; // elements with a finite, positive estimate
    int    n_scanned;      // elements iterated
};

// Scans the elements of theModel's domain.
//   useTangent : getTangentStiff() instead of getInitialStiff().
//   lumping    : lumping of the element mass.
// Under a parallel build the limits are reduced with MPI_MIN over the ranks.
//
// An element without own mass, in a model whose inertia is nodal, yields no
// eigenvalue and contributes no bound.
CTSResult computeCriticalTimeStep(AnalysisModel *theModel,
                                  bool useTangent,
                                  CTSLumping lumping);

// Lumps the consistent mass M of one element into mdiag (sized M.noRows()). M is
// copied before the return, so the caller may fetch the stiffness afterwards: many
// elements return the same shared static matrix from getMass() and getInitialStiff().
void lumpElementMass(Element *ele, const Matrix &M, CTSLumping lumping,
                     double *mdiag);

// Largest generalized eigenvalue of K v = lambda M v for one element with diagonal
// mass mdiag. Returns < 0 on a size mismatch or eigensolver failure.
double elementLambdaMax(Element *ele, bool useTangent,
                        const double *mdiag, int n);

// Sum of the stiffness-proportional Rayleigh coefficients betaK + betaK0 + betaKc,
// each clamped at 0. All three act through a stiffness-like matrix and shrink the
// explicit step identically at the initial state.
double stiffnessRayleighBeta(Element *ele);

// True if every entry of v is finite.
bool vectorIsFinite(const Vector &v);

#endif
