// Ladruno WP-248: `model testNDMaterial` — an nD material tester with explicit commit.
// See TclNDMaterialTester.h and Ladruno_implementation/TclNDMaterialTester_guide.md.
//
//   ndTest material $matTag $view      copy nD material $matTag in $view; returns the order n
//   ndTest order                       the order n of the material under test
//   ndTest strain e1 ... en            set the trial strain (all n components)
//   ndTest mixed $kinds v1 ... vn ?-tol t? ?-maxIter k?
//                                      $kinds is n letters, e (strain given) or s (stress given);
//                                      Newton on the material tangent finds the free strains;
//                                      returns the number of iterations
//   ndTest getStrain | getStress       n numbers
//   ndTest getTangent                  n*n numbers, row by row
//   ndTest commit                      commitState, and remember the committed strain
//   ndTest revert                      revertToLastCommit, then the committed strain as the trial,
//                                      so getStrain/getStress read the committed state
//
// Nothing commits unless `ndTest commit` is called, so trials are free.

#include <TclNDMaterialTester.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

#include <Matrix.h>
#include <NDMaterial.h>
#include <Vector.h>

static TclNDMaterialTester *theTester = 0;
static NDMaterial *theMaterial = 0;
static int theOrder = 0;
static Vector theCommitted;  // the strain at the last ndTest commit (zero at ndTest material)

static int ndTestError(Tcl_Interp *interp, const std::string &message)
{
  std::string text = "ndTest: " + message;
  Tcl_SetResult(interp, const_cast<char *>(text.c_str()), TCL_VOLATILE);
  return TCL_ERROR;
}

static int ndTestNumber(Tcl_Interp *interp, TCL_Char *token, double &value, const char *what)
{
  if (Tcl_GetDouble(interp, token, &value) != TCL_OK || !std::isfinite(value))
    return ndTestError(interp, std::string(what) + " \"" + token + "\" is not a finite number");
  return TCL_OK;
}

static void ndTestAppend(std::string &out, double v)
{
  char buffer[32];
  std::snprintf(buffer, sizeof buffer, "%.17g", v);
  if (!out.empty()) out += ' ';
  out += buffer;
}

static int ndTestVector(Tcl_Interp *interp, const Vector &v)
{
  std::string out;
  for (int i = 0; i < v.Size(); i++) ndTestAppend(out, v(i));
  Tcl_SetResult(interp, const_cast<char *>(out.c_str()), TCL_VOLATILE);
  return TCL_OK;
}

static int ndTestSetStrain(Tcl_Interp *interp, const Vector &e)
{
  if (theMaterial->setTrialStrain(e) < 0)
    return ndTestError(interp, "the material refused the trial strain");
  return TCL_OK;
}

// Mixed control: strain components marked 'e' are imposed; components marked 's'
// are solved so that their stress equals the given value. Newton on K_ss.
static int ndTestMixed(Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  const int n = theOrder;
  if (argc < 3 + n)
    return ndTestError(interp, "want: ndTest mixed kinds v1 ... v" + std::to_string(n) + " ?-tol t? ?-maxIter k?");
  const char *kinds = argv[2];
  if ((int)std::strlen(kinds) != n)
    return ndTestError(interp, std::string("kinds \"") + kinds + "\" must have " + std::to_string(n) + " letters, e or s");
  Vector target(n);
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (kinds[i] != 'e' && kinds[i] != 's')
      return ndTestError(interp, std::string("kinds \"") + kinds + "\": each letter is e (strain) or s (stress)");
    if (kinds[i] == 's') m++;
    if (ndTestNumber(interp, argv[3 + i], target(i), "value") != TCL_OK) return TCL_ERROR;
  }
  double tol = 1.0e-10;
  int maxIter = 25;
  for (int a = 3 + n; a < argc; a += 2) {
    if (a + 1 >= argc) return ndTestError(interp, std::string("option ") + argv[a] + " needs a value");
    if (std::strcmp(argv[a], "-tol") == 0) {
      if (ndTestNumber(interp, argv[a + 1], tol, "-tol") != TCL_OK) return TCL_ERROR;
      if (!(tol > 0.0)) return ndTestError(interp, "-tol must be positive");
    } else if (std::strcmp(argv[a], "-maxIter") == 0) {
      if (Tcl_GetInt(interp, argv[a + 1], &maxIter) != TCL_OK || maxIter < 1 || maxIter > 1000)
        return ndTestError(interp, "-maxIter must be a whole number from 1 to 1000");
    } else {
      return ndTestError(interp, std::string("unknown option ") + argv[a]);
    }
  }

  // start from the current trial strain, with the imposed components set
  Vector e(theMaterial->getStrain());
  if (e.Size() != n) return ndTestError(interp, "the material's strain has the wrong size");
  for (int i = 0; i < n; i++)
    if (kinds[i] == 'e') e(i) = target(i);

  Matrix Kss(m > 0 ? m : 1, m > 0 ? m : 1);
  Vector r(m > 0 ? m : 1), d(m > 0 ? m : 1);
  double rmax = 0.0;
  for (int iter = 0; iter <= maxIter; iter++) {
    if (ndTestSetStrain(interp, e) != TCL_OK) return TCL_ERROR;
    if (m == 0) {
      Tcl_SetObjResult(interp, Tcl_NewIntObj(0));
      return TCL_OK;
    }
    const Vector &sigma = theMaterial->getStress();
    double scale = 0.0;
    rmax = 0.0;
    for (int i = 0, k = 0; i < n; i++) {
      scale = std::fmax(scale, std::fabs(sigma(i)));
      if (kinds[i] != 's') continue;
      scale = std::fmax(scale, std::fabs(target(i)));
      r(k) = sigma(i) - target(i);
      rmax = std::fmax(rmax, std::fabs(r(k)));
      k++;
    }
    if (!std::isfinite(rmax)) return ndTestError(interp, "the stress is not finite");
    if (rmax <= tol * scale) {
      Tcl_SetObjResult(interp, Tcl_NewIntObj(iter));
      return TCL_OK;
    }
    if (iter == maxIter) break;
    const Matrix &K = theMaterial->getTangent();
    for (int i = 0, a = 0; i < n; i++) {
      if (kinds[i] != 's') continue;
      for (int j = 0, b = 0; j < n; j++) {
        if (kinds[j] != 's') continue;
        Kss(a, b++) = K(i, j);
      }
      a++;
    }
    if (Kss.Solve(r, d) < 0)
      return ndTestError(interp, "the tangent of the stress-controlled components is singular");
    for (int i = 0, k = 0; i < n; i++)
      if (kinds[i] == 's') e(i) -= d(k++);
  }
  char buffer[160];
  std::snprintf(buffer, sizeof buffer, "mixed control did not converge in %d iterations (stress residual %.6g)", maxIter, rmax);
  return ndTestError(interp, buffer);
}

static int ndTestCommand(ClientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theTester == 0) return ndTestError(interp, "the tester has been removed");
  if (argc < 2)
    return ndTestError(interp, "want: ndTest material|order|strain|mixed|getStrain|getStress|getTangent|commit|revert ...");
  const char *sub = argv[1];

  if (std::strcmp(sub, "material") == 0) {
    if (argc != 4) return ndTestError(interp, "want: ndTest material matTag view");
    int tag;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) return ndTestError(interp, std::string("matTag \"") + argv[2] + "\" is not an integer");
    NDMaterial *orig = OPS_getNDMaterial(tag);
    if (orig == 0) return ndTestError(interp, "no nDMaterial with tag " + std::to_string(tag));
    NDMaterial *copy = orig->getCopy(argv[3]);
    if (copy == 0) return ndTestError(interp, std::string("nDMaterial ") + std::to_string(tag) + " (" + orig->getType() + ") has no " + argv[3] + " view");
    const int n = copy->getStrain().Size();
    if (n < 1) {
      delete copy;
      return ndTestError(interp, std::string("the ") + argv[3] + " view of nDMaterial " + std::to_string(tag) + " has no strain components");
    }
    if (theMaterial != 0) delete theMaterial;
    theMaterial = copy;
    theOrder = n;
    theCommitted = Vector(n);
    Tcl_SetObjResult(interp, Tcl_NewIntObj(n));
    return TCL_OK;
  }

  if (theMaterial == 0) return ndTestError(interp, "no material under test: call ndTest material matTag view first");

  if (std::strcmp(sub, "order") == 0) {
    Tcl_SetObjResult(interp, Tcl_NewIntObj(theOrder));
    return TCL_OK;
  }
  if (std::strcmp(sub, "strain") == 0) {
    if (argc != 2 + theOrder) return ndTestError(interp, "want: ndTest strain e1 ... e" + std::to_string(theOrder));
    Vector e(theOrder);
    for (int i = 0; i < theOrder; i++)
      if (ndTestNumber(interp, argv[2 + i], e(i), "strain") != TCL_OK) return TCL_ERROR;
    return ndTestSetStrain(interp, e);
  }
  if (std::strcmp(sub, "mixed") == 0) return ndTestMixed(interp, argc, argv);
  if (std::strcmp(sub, "getStrain") == 0) return ndTestVector(interp, theMaterial->getStrain());
  if (std::strcmp(sub, "getStress") == 0) return ndTestVector(interp, theMaterial->getStress());
  if (std::strcmp(sub, "getTangent") == 0) {
    const Matrix &K = theMaterial->getTangent();
    std::string out;
    for (int i = 0; i < K.noRows(); i++)
      for (int j = 0; j < K.noCols(); j++) ndTestAppend(out, K(i, j));
    Tcl_SetResult(interp, const_cast<char *>(out.c_str()), TCL_VOLATILE);
    return TCL_OK;
  }
  if (std::strcmp(sub, "commit") == 0) {
    if (theMaterial->commitState() < 0) return ndTestError(interp, "the material refused to commit");
    theCommitted = theMaterial->getStrain();
    return TCL_OK;
  }
  if (std::strcmp(sub, "revert") == 0) {
    if (theMaterial->revertToLastCommit() < 0) return ndTestError(interp, "the material refused to revert");
    return ndTestSetStrain(interp, theCommitted);
  }
  return ndTestError(interp, std::string("unknown subcommand ") + sub);
}

TclNDMaterialTester::TclNDMaterialTester(Domain &theDomain, Tcl_Interp *interp)
  : TclModelBuilder(theDomain, interp, 1, 1), theInterp(interp)
{
  Tcl_CreateCommand(interp, "ndTest", ndTestCommand, (ClientData)NULL, NULL);
  theTester = this;
}

TclNDMaterialTester::~TclNDMaterialTester()
{
  theTester = 0;
  if (theMaterial != 0) {
    delete theMaterial;
    theMaterial = 0;
  }
  theOrder = 0;
  Tcl_DeleteCommand(theInterp, "ndTest");
}
