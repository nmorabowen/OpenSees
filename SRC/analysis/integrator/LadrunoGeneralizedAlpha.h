/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno fork — N. Mora-Bowen
// ADR-52 W3-I2: sensitivity-carrying (DDM) subclass of GeneralizedAlpha.
//
// Strict superset of LadrunoHHT: the Chung-Hulbert generalized-alpha method
// carries TWO spectral parameters — alphaF weights the stiffness/damping
// (evaluated at Ualpha / Ualphadot) and alphaM weights the inertia (evaluated at
// Ualphadotdot) — where HHT has a single alpha on K/C and full-step M. The DDM
// seam therefore generalizes LadrunoHHT's:
//   - the damping multiplicator picks up the alphaF-weighting of d(Ualphadot)/dh,
//   - the mass multiplicator picks up the alphaM-weighting of d(Ualphadotdot)/dh,
//   - dC/dh acts on Ualphadot, dM/dh on Ualphadotdot,
//   - the extra element-only stiffness term is -K*(1-alphaF)*dU_n (the (1-alphaF)
//     part of d(Ualpha)/dh).
// At alphaM == alphaF == 1 every weighting collapses and this reduces to Newmark
// DDM (average acceleration). Because U/Udot/Udotdot evolve by the SAME Newmark
// relations, saveSensitivity / commitSensitivity / formSensitivityRHS /
// formIndependentSensitivityRHS / computeSensitivities are identical to Newmark.
//
// See Ladruno_implementation/52_ladruno_integrator_strengthening_adr.md (W3-I2).

#ifndef LadrunoGeneralizedAlpha_h
#define LadrunoGeneralizedAlpha_h

#include <GeneralizedAlpha.h>
#include <Vector.h>

class DOF_Group;
class FE_Element;

class LadrunoGeneralizedAlpha : public GeneralizedAlpha
{
public:
    // constructors (same signatures as GeneralizedAlpha; own classTag)
    LadrunoGeneralizedAlpha();
    LadrunoGeneralizedAlpha(double alphaM, double alphaF);
    LadrunoGeneralizedAlpha(double alphaM, double alphaF, double beta, double gamma);

    // destructor
    ~LadrunoGeneralizedAlpha();

    // residual hooks branch on sensitivityFlag (DDM); the non-sensitivity path
    // delegates to the base so a plain transient run is byte-identical.
    int formEleResidual(FE_Element *theEle);
    int formNodUnbalance(DOF_Group *theDof);

    // tangent hooks branch on sensTangentFlag: the non-sensitivity path delegates
    // to the base (byte-identical primal). The sensitivity path emits the
    // PRIMAL-CONSISTENT Jacobian (M-coef c3, NOT the base's alphaM*c3) so the DDM
    // solve uses the true ∂R/∂U of the as-implemented generalized-alpha residual
    // (inertia at Udotdot). See the .cpp header for why the base tangent is
    // inconsistent and why the sensitivity solve must re-form its own tangent.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);

    const char *getClassType(void) const;
    void Print(OPS_Stream &s, int flag = 0);

    // AddingSensitivity:BEGIN //////////////////////////////////
    int revertToStart();
    int formSensitivityRHS(int gradNum);
    int formIndependentSensitivityRHS();
    int saveSensitivity(const Vector &v, int gradNum, int numGrads);
    int commitSensitivity(int gradNum, int numGrads);
    int computeSensitivities();
    // AddingSensitivity:END ////////////////////////////////////

protected:

private:
    // DDM workspace — same roles as Newmark's protected sensitivity block.
    int sensitivityFlag;
    int gradNumber;
    Vector *massMatrixMultiplicator;
    Vector *dampingMatrixMultiplicator;
    int assemblyFlag;
    Vector independentRHS;
    int sensTangentFlag;   // 1 while re-forming the consistent tangent for the DDM solve
};

#endif
