/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ladruno (ADR-87 D2) — mutation gates: making "the tests are green" falsifiable.
//
// WHAT
//   A mutation build deliberately sabotages ONE family's physics. CI then runs
//   that family's test suite and requires it to turn RED. A suite that stays
//   green while the physics is deleted is not testing the physics — that is a
//   defect of the same severity as a wrong stress, and this header is the
//   machinery that proves it either way.
//
// WHY THIS SHAPE (three deliberate choices)
//
//   1. ALWAYS COMPILED, never #ifdef'd out. The mutation call sites are real
//      C++ in every build; only the *mode constant* changes. So a refactor
//      that renames the vector a gate reads breaks the build immediately,
//      instead of silently disabling the gate and leaving a green CI lane that
//      proves nothing. Scaffolding that can rot into a no-op is worse than no
//      scaffolding: it reports safety it is not delivering.
//
//   2. ZERO COST WHEN OFF. With mode == LADRUNO_MUT_NONE every branch below is
//      `if (0 == k)`, which every optimizer folds away. The default build
//      carries no mutation code in its object files, so this cannot change a
//      shipped answer. (Verify, don't assume: `tests/test_mutation_framework.py`
//      pins that a default build reports no active mutation at runtime.)
//
//   3. MUTATE THE ELEMENT<->SOLVER CONTRACT, NOT THE INTERNALS. Gates are
//      placed in the Element API accessors (getTangentStiff, getInitialStiff,
//      getResistingForce, getResistingForceIncInertia), which is everything an
//      analysis can observe about an element. That makes a gate
//      formulation-agnostic — LadrunoBrick's six assembly paths (std/bbar,
//      URI, physical-hourglass, SSP, EAS, finite, hypo) are covered by four
//      call sites, and a NEW formulation is covered the day it is added,
//      because it cannot reach the solver except through these functions.
//
// MODES
//   ZERO   internal force -> 0. The canonical "replace the physics with
//          return 0". Any test that checks equilibrium, a reaction, or a
//          displacement must fail.
//   SCALE  internal force *= 1.5. The DIAGNOSTIC mode: it yields a
//          plausible-but-wrong answer rather than an obviously broken one, so
//          it separates suites that check values from suites that only check
//          "it ran without throwing" — precisely the failure mode ADR-87 was
//          written about. A suite green under SCALE is a smoke test.
//   IDENT  tangent -> identity, residual untouched. Newton usually still
//          converges to the RIGHT answer (slowly), so most tests correctly
//          stay green; this mode exists to ask whether anything at all pins
//          the tangent (iteration counts, eigenvalues, convergence rate).
//          Read a low IDENT score as "no tangent coverage", not as a bug.
//
// USAGE (one line per accessor; FAMILY is a bare token, e.g. CONTINUUM)
//   const Vector &Foo::getResistingForce(void) {
//     formResidAndTangent(0);
//     LADRUNO_MUTATE_FORCE(CONTINUUM, resid);   // ADR-87 D2 gate
//     ...
//   }
//
// ADDING A FAMILY: add its LADRUNO_MUTATE_<FAMILY> default below, put the two
// macros in its Element API accessors, register it in
// Ladruno_scripts/mutation_gate.py, and record the score in the family's
// verification manifest (ADR-87 D1). All four steps, or the family is not
// `shipped`.
//
// Grep the live gate list with:  grep -rn "LADRUNO_MUTATE_" SRC/

#ifndef Ladruno_mutation_h
#define Ladruno_mutation_h

#include <Vector.h>
#include <Matrix.h>

// ---- modes ---------------------------------------------------------------
#define LADRUNO_MUT_NONE  0
#define LADRUNO_MUT_ZERO  1
#define LADRUNO_MUT_SCALE 2
#define LADRUNO_MUT_IDENT 3

// ---- per-family selector --------------------------------------------------
// Set by the build system, e.g. -DLADRUNO_MUTATE_CONTINUUM=1. Anything not
// set by the build defaults to NONE, so a default build is unmutated.
#ifndef LADRUNO_MUTATE_CONTINUUM
#  define LADRUNO_MUTATE_CONTINUUM LADRUNO_MUT_NONE
#endif
#ifndef LADRUNO_MUTATE_UP
#  define LADRUNO_MUTATE_UP LADRUNO_MUT_NONE
#endif
#ifndef LADRUNO_MUTATE_CONTACT
#  define LADRUNO_MUTATE_CONTACT LADRUNO_MUT_NONE
#endif
#ifndef LADRUNO_MUTATE_EXPLICIT
#  define LADRUNO_MUTATE_EXPLICIT LADRUNO_MUT_NONE
#endif
#ifndef LADRUNO_MUTATE_SANISAND
#  define LADRUNO_MUTATE_SANISAND LADRUNO_MUT_NONE
#endif
// ADR 91: the shell section stiffness-modifier decorator. Gated at the SECTION
// resultant/tangent, not at an element -- the whole feature IS a constitutive
// transform, so the section accessors are its only physics surface.
#ifndef LADRUNO_MUTATE_SHELLMOD
#  define LADRUNO_MUTATE_SHELLMOD LADRUNO_MUT_NONE
#endif

namespace LadrunoMutation {

// The SCALE factor is deliberately NOT near 1.0: a 50% error must be caught by
// any test that checks a value to engineering tolerance, so a green result
// under SCALE means the value is unchecked, never "the tolerance was tight".
const double SCALE_FACTOR = 1.5;

// True iff any family is mutated. Reported by the `ladrunoMutation` verb so a
// harness can refuse to interpret results from a mutant binary as physics.
inline bool anyActive(void)
{
  return (LADRUNO_MUTATE_CONTINUUM != LADRUNO_MUT_NONE) ||
         (LADRUNO_MUTATE_UP        != LADRUNO_MUT_NONE) ||
         (LADRUNO_MUTATE_CONTACT   != LADRUNO_MUT_NONE) ||
         (LADRUNO_MUTATE_EXPLICIT  != LADRUNO_MUT_NONE) ||
         (LADRUNO_MUTATE_SANISAND  != LADRUNO_MUT_NONE) ||
         (LADRUNO_MUTATE_SHELLMOD  != LADRUNO_MUT_NONE);
}

// Sabotage an internal-force vector. Placed AFTER the force is formed and
// BEFORE element loads/inertia are folded in, so the mutation hits the
// constitutive path only — an inertia- or load-driven test that survives is
// telling the truth about what it covers.
inline void force(int mode, Vector &v)
{
  if (mode == LADRUNO_MUT_ZERO)
    v.Zero();
  else if (mode == LADRUNO_MUT_SCALE)
    v *= SCALE_FACTOR;
  // IDENT is a tangent-only mode: the residual is deliberately left exact so
  // the run still converges and the gate measures tangent coverage alone.
}

// Sabotage a stiffness matrix. IDENT replaces K with the identity, which keeps
// the matrix non-singular (the analysis proceeds and Newton can still find the
// right root) so the gate measures whether anything CHECKS the tangent rather
// than merely whether the solve survives.
inline void tangent(int mode, Matrix &m)
{
  if (mode == LADRUNO_MUT_ZERO)
    m.Zero();
  else if (mode == LADRUNO_MUT_SCALE)
    m *= SCALE_FACTOR;
  else if (mode == LADRUNO_MUT_IDENT) {
    m.Zero();
    const int n = (m.noRows() < m.noCols()) ? m.noRows() : m.noCols();
    for (int i = 0; i < n; i++)
      m(i, i) = 1.0;
  }
}

} // namespace LadrunoMutation

// ---- call-site macros -----------------------------------------------------
// FAMILY is a bare token (CONTINUUM, UP, ...). Token-pasting resolves it to
// that family's mode constant, so an unknown family is a COMPILE error rather
// than a silently inactive gate.
#define LADRUNO_MUTATE_FORCE(FAMILY, VEC) \
  LadrunoMutation::force(LADRUNO_MUTATE_##FAMILY, (VEC))

#define LADRUNO_MUTATE_TANGENT(FAMILY, MAT) \
  LadrunoMutation::tangent(LADRUNO_MUTATE_##FAMILY, (MAT))

#endif // Ladruno_mutation_h
