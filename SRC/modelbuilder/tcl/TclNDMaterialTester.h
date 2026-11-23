// Ladruno WP-248: `model testNDMaterial` — an nD material tester with explicit commit.
//
// The stock testers (`testUniaxial`, `testPlaneStress`) commit on every N-th strain
// call and cannot revert, so a stress-controlled step has to be searched with a fixed
// number of trials. This tester holds a copy of one NDMaterial in a requested view
// (ThreeDimensional, PlaneStrain, PlaneStress, ...), takes trial strains or a mixed
// strain/stress control solved by Newton on the material tangent, and commits or
// reverts only when told to. It registers one Tcl command, `ndTest`; see
// Ladruno_implementation/TclNDMaterialTester_guide.md.

#ifndef TclNDMaterialTester_h
#define TclNDMaterialTester_h

#include <TclModelBuilder.h>
#include <tcl.h>

class TclNDMaterialTester : public TclModelBuilder
{
  public:
    TclNDMaterialTester(Domain &theDomain, Tcl_Interp *interp);
    ~TclNDMaterialTester();

  private:
    Tcl_Interp *theInterp;
};

#endif
