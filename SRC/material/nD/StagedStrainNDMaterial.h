/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
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

// Authors: Nicolas Mora Bowen, Guppi (Ladruño)
// Created: 06/2026
//
// Description: StagedStrainNDMaterial — the SMALL-STRAIN, dimension-general,
// AUTO-CAPTURING staged-activation wrapper. It is the small-strain sibling of
// InitDefGrad (StagedDefGrad): it makes a small-strain continuum element BORN
// STRESS-FREE at the current deformed geometry when appended mid-stage (staged
// construction). User-facing command: `nDMaterial StagedStrain`.
//
// WHY A NEW CLASS (and not the upstream InitStrain). InitStrainNDMaterial (Petracca,
// 2024) is a fine FIXED, ADDITIVE prestrain (thermal / lack-of-fit), but for staged
// construction it bites three ways:
//   1. 3D ONLY — getOrder() is hardcoded to 6 and getCopy("ThreeDimensional") is the
//      only path, so a 2D / PlaneStrain / PlaneStress / AxiSymmetric element cannot
//      use it; worse, a non-3D inner trips its exit(-1) and KILLS the interpreter.
//   2. NO AUTO-CAPTURE — ε0 is fixed at construction, so one tag cannot birth a field
//      of elements with different birth strains.
//   3. Its intent is to ADD a known prestrain, not to SUBTRACT the captured birth
//      strain for stress-free birth.
// StagedStrain fixes all three: dimension-general, auto-capturing, graceful, and
// subtractive (re-references to birth). InitStrain is left untouched for prestrain.
//
// THE MECHANIC (additive re-reference; small strain ⇒ exact, no objectivity needed).
// At activation the wrapper captures the birth strain
//
//     ε0 = ε at the first setTrialStrain   ( = B·u0, the committed strain at birth )
//
// and thereafter feeds the inner material the RELATIVE strain
//
//     ε_rel = ε − ε0    ( deformation measured from the stress-free birth state )
//
// At birth ε_rel = 0 ⇒ the inner (started pristine) returns ZERO stress and ZERO
// plastic history — the element is born GENUINELY VIRGIN, not "stress zeroed but
// history carried". It then accrues stress only from deformation BEYOND birth. The
// tangent is an EXACT passthrough (∂σ/∂ε = ∂σ/∂(ε−ε0) since ε0 is constant), so no
// finite-difference tangent gate is needed (unlike the finite InitDefGrad). This is
// EXACT for small strain; for finite strain the additive offset is non-objective —
// use StagedDefGrad (InitDefGrad) there.
//
// DIMENSION GENERALITY (the thing InitStrain skipped). The wrapper is agnostic to the
// material "type": ε0 is sized to inner->getOrder() (3 for PlaneStrain/PlaneStress,
// 4 for AxiSymmetric, 6 for ThreeDimensional, …) and getCopy(type) adapts the inner
// to whatever dimensional view the element requests, so ONE tag serves 2D and 3D
// elements alike — each copy auto-captures at its own order.
//
// AUTO-CAPTURE TRIGGER. The wrapper material instance is created fresh per element
// (getCopy per Gauss point), so its first setTrialStrain carries the committed birth
// strain (addElement→setDomain→update evaluates at the committed deformed state).
// Each appended element captures its own ε0 ⇒ multi-lift falls out for free, no
// InitialStateAnalysis-style global toggle. Degenerates to identity at t=0 (ε0=0).
//
// COMPOSABILITY (staged birth + prestrain). StagedStrain and the upstream InitStrain
// compose by nesting, with StagedStrain OUTERMOST (element-facing):
//
//     element → StagedStrain( InitStrain( realMaterial, ε_pre ) )
//
// At birth StagedStrain feeds InitStrain ε_rel=0, InitStrain adds ε_pre, so the real
// material sees ε_pre ⇒ σ(ε_pre): born carrying EXACTLY the prestress at the deformed
// geometry, no inherited geometric stress. (Nest the other way and the prestrain is
// captured into ε0 and cancels — so StagedStrain must be the outer wrapper.)
//
// OPTIONS:  -noInit            pass-through (ε0 ≡ 0; bit-identical to the bare inner).
//           -eps0 e1 e2 …      explicit birth strain (skip auto-capture); component
//                              count must match the inner order (subtract convention,
//                              i.e. the OPPOSITE sign of InitStrain's additive prestrain).
//
// See Ladruno_implementation/staged_deformation_gradiend.md (the staged-activation
// design note) and constraints_reference_position.md (the SP/MP/element reference-
// frame trichotomy).

#ifndef StagedStrainNDMaterial_h
#define StagedStrainNDMaterial_h

#include <NDMaterial.h>
#include <Vector.h>

class StagedStrainNDMaterial : public NDMaterial
{
 public:
  // active=false  -> pure pass-through (ε0 ≡ 0; the -noInit opt-out).
  // eps0Given!=0  -> explicit birth strain (size must equal the inner order); skips
  //                  the auto-capture. Otherwise ε0 is captured on the FIRST setTrialStrain.
  StagedStrainNDMaterial(int tag, NDMaterial &inner, bool active = true,
                         const Vector *eps0Given = 0);
  StagedStrainNDMaterial();
  ~StagedStrainNDMaterial();

  const char *getClassType(void) const { return "StagedStrainNDMaterial"; }

  // --- strain seam (any dimensional order) ---
  int setTrialStrain(const Vector &v);
  int setTrialStrain(const Vector &v, const Vector &rate);
  int setTrialStrainIncr(const Vector &v);
  int setTrialStrainIncr(const Vector &v, const Vector &rate);

  const Vector &getStress(void);
  const Matrix &getTangent(void);
  const Matrix &getInitialTangent(void);
  const Vector &getStrain(void);      // RELATIVE strain (ε − ε0), the physical state
  double getRho(void);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  NDMaterial *getCopy(void);
  NDMaterial *getCopy(const char *type);   // adapts the inner to the requested view
  const char *getType(void) const;         // delegates to the inner
  int getOrder(void) const;                // delegates to the inner

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  Response *setResponse(const char **argv, int argc, OPS_Stream &output);
  int getResponse(int responseID, Information &matInfo);

  // sensitivity (DDM) — parity with InitStrainNDMaterial
  const Vector &getStressSensitivity(int gradIndex, bool conditional);
  int commitSensitivity(const Vector &depsdh, int gradIndex, int numGrads);

 private:
  NDMaterial *theMaterial;   // inner small-strain material (deep copy)

  bool   active;             // false ⇒ pass-through (no re-referencing)
  bool   captured;           // ε0 has been snapshot
  bool   eps0Explicit;       // ε0 supplied at construction (revertToStart keeps it)
  Vector eps0;               // birth strain (sized to inner order)

  Vector relStrain;          // scratch: ε − ε0 forwarded to the inner
  Vector totalStrain;        // last total strain seen (for the totalStrain response)

  void sizeBuffers(int n);   // size ε0/relStrain/totalStrain to a given order
};

#endif
