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
#define LADRUNO_MATERIAL_REFUSED (-33086)

#endif
