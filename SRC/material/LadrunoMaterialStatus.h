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

// Ladruno (ADR-86b): the ONE material return code a fork element propagates.
//
// WHY THIS IS NOT JUST `return -1`.
//
// OpenSees has no convention distinguishing the two things a material can mean
// by a non-zero `setTrialStrain`:
//
//   (a) "an inner iteration missed, but here is my best state" -- transient,
//       recoverable, and the GLOBAL Newton is the right arbiter. The fork
//       already has a hard-won rule for this case (LEDGER_quirks, ADR-33/34):
//       return 0 with the last iterate plus a loud warning, because returning
//       a failure code makes softening analyses FRAGILE -- a fixed-increment
//       run just dies at a kink the global Newton would have walked through.
//       `ASDConcrete3DMaterial::setTrialStrain` returns negative in exactly
//       this sense, and MEASURED (WP-86b): making `LadrunoBrick::update()`
//       propagate any `< 0` killed `test_ladrunoBrick_asdconcrete_bend.py`'s
//       two mesh-objectivity gates at load factor 605, on a run that had been
//       green for months.
//
//   (b) "I did NOT integrate this increment; my committed state is unchanged;
//       cut the step." That is a different statement, and it is the one the
//       ADR-86b substep cap makes. Nothing recoverable happened, there is no
//       best-effort state to accept, and swallowing it reproduces exactly the
//       silence ADR-90 GATE U ran into (0 of 80 subdivisions used across six
//       legs while single steps ran for 34 minutes).
//
// So (b) gets its own value and elements propagate ONLY that. Every other
// non-zero code keeps whatever treatment it already had, which is why adopting
// this cost zero behaviour change anywhere else in the tree.
//
// The value is negative (so `< 0` tests still see a failure) and far from the
// small integers materials return by hand, so it cannot be produced by
// accident. It is NOT a class tag and does not belong in classTags.h.
//
// See Ladruno_implementation/86_ladruno_sanisand_handoff.md and LEDGER_quirks.

#ifndef LadrunoMaterialStatus_h
#define LadrunoMaterialStatus_h

// Returned by a Ladruno material whose update did not integrate the increment.
// The committed state is guaranteed untouched; the caller must fail the step.
//
// `constexpr`, not a macro: typed, scoped, and it cannot collide with anything a
// later header spells the same way. The VALUE is arbitrary -- it only has to be
// negative (so pre-existing `< 0` tests still see a failure) and far from the
// small integers materials return by hand. It is NOT a class tag, is not
// registered in classTags.h, and nothing may derive one from it.
constexpr int LADRUNO_MATERIAL_REFUSED = -33086;

// ==========================================================================
//  Ladruno WP-99 (F7): the COMMIT-TIME refusal seam.
//
//  The sentinel above only works at the TRIAL, and only on an element that
//  forwards `setTrialStrain`'s return code. At COMMIT there is no such path at
//  all: `Domain::commit()` is `elePtr->commitState();` with the return value
//  dropped, for EVERY element, fork or vanilla. A material that discovers at
//  commitState() that it cannot integrate the step therefore had nowhere to
//  say so -- measured on the TIMs strip as 25.9 M capped commits and a
//  straight-line load-settlement curve, every step reported converged.
//
//  So the refusal goes around the element instead of through it. A material
//  calls ladrunoNoteCommitRefusal() from its commitState(); Domain::commit()
//  checks the count after its element loop and aborts the commit with a
//  negative return, which AnalysisModel::commitDomain() turns into -2
//  (AnalysisModel.cpp:656-659) and every analysis class turns into a failed
//  step: StaticAnalysis.cpp:213-222 and DirectIntegrationAnalysis.cpp:259-267
//  both revert and `return -4`, and
//  VariableTimeStepDirectIntegrationAnalysis.cpp:137-140 sets result = -4 and
//  falls into its revert-and-subdivide branch. That is element-independent by
//  construction, which is the whole point: a discarding element (Brick,
//  BbarBrick, SSPquad, ...) cannot swallow it.
//
//  WHY NOT SIMPLY PROPAGATE `elePtr->commitState()`'s RETURN? Because ADR-33/34
//  forbids it: ASDConcrete3D and friends return negative "best-state" codes
//  from a commit that is perfectly valid, and failing the step on those was
//  MEASURED to break mesh-objectivity gates that had been green for months.
//  The rule the fork settled on is that only a DECLARED refusal fails a step,
//  never any nonzero code -- and at commit there is no sentinel-filtering
//  element in the path to apply that rule, so the declaration has to arrive
//  out of band. This counter is that declaration.
//
//  Header-only, with a function-local static in an inline function (one
//  instance across all TUs by the ODR), so no library-dependency edge is added
//  between SRC/domain and SRC/material.
//
//  NOT thread-safe, deliberately: it is written from the commit phase, which is
//  serial in every analysis class today. If ADR-75b Lane 3 ever threads
//  commitState(), this becomes an atomic.
// ==========================================================================
inline int &ladrunoCommitRefusalCounter(void)
{
    static int nCommitRefusals = 0;
    return nCommitRefusals;
}

// Called by a material whose commitState() could not integrate the step.
inline void ladrunoNoteCommitRefusal(void) { ++ladrunoCommitRefusalCounter(); }

// How many integration points refused the commit now being assembled.
inline int ladrunoPendingCommitRefusals(void) { return ladrunoCommitRefusalCounter(); }

// Called by Domain::commit() once it has acted on them.
inline void ladrunoClearCommitRefusals(void) { ladrunoCommitRefusalCounter() = 0; }

#endif
