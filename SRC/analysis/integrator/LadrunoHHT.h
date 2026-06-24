/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno fork — N. Mora-Bowen
// ADR-52 W3-I2: sensitivity-carrying (DDM) subclass of HHT.
//
// LadrunoHHT reuses the entire HHT algorithm (newStep / update / commit /
// tangent — inherited unchanged) and adds the Direct Differentiation Method
// (DDM) sensitivity seam that, in vanilla OpenSees, exists only on Newmark /
// DisplacementControl / ArcLength / MinUnbalDispNorm. This unblocks
// reliability / fragility / FORM on the numerically-damped HHT integrator.
//
// The seam mirrors Newmark's five-method pattern, BUT the dynamic residual is
// evaluated at the HHT alpha-intermediate state (Ualpha / Ualphadot) instead of
// the end-of-step state. Because HHT's U / Udot / Udotdot evolve by the SAME
// Newmark relations, saveSensitivity / commitSensitivity / formSensitivityRHS /
// formIndependentSensitivityRHS / computeSensitivities are identical to Newmark.
// Only formEleResidual / formNodUnbalance change:
//   - the damping multiplicator picks up the alpha-weighting of d(Ualphadot)/dh,
//   - the C-matrix sensitivity acts on Ualphadot (not Udot),
//   - an EXTRA element-only stiffness term  -K*(1-alpha)*dU_n  appears (the
//     (1-alpha) part of d(Ualpha)/dh); it vanishes for Newmark (alpha == 1).
// The mass terms are unchanged (M acts at the full step in HHT).
//
// See Ladruno_implementation/52_ladruno_integrator_strengthening_adr.md (W3-I2).

#ifndef LadrunoHHT_h
#define LadrunoHHT_h

#include <HHT.h>
#include <Vector.h>

class DOF_Group;
class FE_Element;

class LadrunoHHT : public HHT
{
public:
    // constructors (same signatures as HHT; register under INTEGRATOR_TAGS_LadrunoHHT)
    LadrunoHHT();
    LadrunoHHT(double alpha);
    LadrunoHHT(double alpha, double beta, double gamma);

    // destructor
    ~LadrunoHHT();

    // residual hooks branch on sensitivityFlag (DDM); the non-sensitivity path
    // delegates to the base so a plain transient run is byte-identical to HHT.
    int formEleResidual(FE_Element *theEle);
    int formNodUnbalance(DOF_Group *theDof);

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
};

#endif
