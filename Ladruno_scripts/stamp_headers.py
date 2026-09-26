#!/usr/bin/env python3
"""Stamp the LADRUNO author/banner header onto every fork-authored source file.

Source of truth for the ASCII art is `Ladruno_scripts/banner_ASCII.txt` (the same
file that drives the runtime splash via patch_banner.py), so the art never drifts
between the console banner and the file headers.

The block is delimited by `// LADRUNO-HEADER-START` / `// LADRUNO-HEADER-END`
markers and is inserted *after* the leading OpenSees/PEER `/* ... */` comment so
upstream attribution stays first. Running again is idempotent: an existing block
is replaced in place (edit the art or credit here, re-run, done). Original line
endings (LF/CRLF) and any BOM are preserved so the stamp produces a clean diff.

    python Ladruno_scripts/stamp_headers.py            # stamp all authored files
    python Ladruno_scripts/stamp_headers.py --check     # report only, exit 1 if any stale
                                                        # or any Ladruno-named SRC source
                                                        # is missing from GLOBS (CI gate)
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# --- credit line (the four authors) ----------------------------------------
CREDIT = [
    "Ladruno — a research fork of OpenSees",
    "Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  "
    "José Abell  ·  Guppi",
]

START = "// LADRUNO-HEADER-START"
END = "// LADRUNO-HEADER-END"
RULE = "// " + "=" * 74

# --- the set of files WE authored (new fork code, not vanilla edits) --------
# Whole directories that are entirely ours, plus specific files in shared dirs.
GLOBS = [
    "SRC/element/ladrunoBrick/*.cpp", "SRC/element/ladrunoBrick/*.h",
    "SRC/element/ladrunoSolidShell/*.cpp", "SRC/element/ladrunoSolidShell/*.h",
    "SRC/element/ladrunoIMKBeam/*.cpp", "SRC/element/ladrunoIMKBeam/*.h",
    "SRC/element/ladrunoEmbeddedRebar/*.cpp", "SRC/element/ladrunoEmbeddedRebar/*.h",
    "SRC/element/ladrunoEmbeddedNode/*.cpp", "SRC/element/ladrunoEmbeddedNode/*.h",
    "SRC/element/ladrunoPlane/*.cpp", "SRC/element/ladrunoPlane/*.h",
    "SRC/element/ladrunoUP/*.cpp", "SRC/element/ladrunoUP/*.h",
    "SRC/element/ladrunoDistributingCoupling/*.cpp", "SRC/element/ladrunoDistributingCoupling/*.h",
    "SRC/element/ladrunoKinematicCoupling/*.cpp", "SRC/element/ladrunoKinematicCoupling/*.h",
    "SRC/element/ladrunoRigidBody/*.cpp", "SRC/element/ladrunoRigidBody/*.h",
    "SRC/element/ladrunoDispBeamColumn/*.cpp", "SRC/element/ladrunoDispBeamColumn/*.h",   # WP-116: was unstamped, so the quirk lint skipped it
    "SRC/element/bezierTriangle/*.cpp", "SRC/element/bezierTriangle/*.h",
    "SRC/element/bezierTetrahedron/*.cpp", "SRC/element/bezierTetrahedron/*.h",
    "SRC/element/solidTransformation/*.cpp", "SRC/element/solidTransformation/*.h",
    "SRC/utility/profiler/*.cpp", "SRC/utility/profiler/*.h",
    "SRC/analysis/analysis/LadrunoComplexEigen.*",
    "SRC/analysis/analysis/LadrunoDampingAssembler.*",
    "SRC/analysis/analysis/LadrunoModalResponse.*",
    "SRC/analysis/analysis/LadrunoModalCombination.*",
    "SRC/system_of_eqn/eigenSOE/FeastEigenSOE.*",
    "SRC/system_of_eqn/eigenSOE/FeastEigenSolver.*",
    "SRC/system_of_eqn/eigenSOE/LadrunoBlockZKernel.*",
    "SRC/system_of_eqn/eigenSOE/LadrunoFeastInnerSolve.*",
    "SRC/system_of_eqn/eigenSOE/LadrunoDistBlockZKernel.*",
    "SRC/system_of_eqn/ladrunoCMS/*.cpp", "SRC/system_of_eqn/ladrunoCMS/*.h",
    "SRC/material/nD/LadrunoJ2.*", "SRC/material/nD/LadrunoJ2Kernel.h",
    "SRC/material/nD/LadrunoJ2Finite.*",
    "SRC/material/nD/LadrunoRCConcrete.*", "SRC/material/nD/LadrunoRCKernel.h",
    "SRC/material/nD/LadrunoRCFiniteStrain.*",
    "SRC/material/LadrunoMaterialStatus.h",
    "SRC/material/nD/LadrunoSANISAND.*",
    "SRC/material/nD/LadrunoSANISAND3D.*",
    "SRC/material/nD/LadrunoSANISANDPlaneStrain.*",
    "SRC/material/nD/LadrunoConcrete3D.*", "SRC/material/nD/LadrunoConcrete3DKernel.h",
    "SRC/material/nD/LadrunoHardening.h",
    "SRC/material/uniaxial/LadrunoUniaxialJ2.*",
    "SRC/material/uniaxial/LadrunoRebarBuckling.*",
    "SRC/material/uniaxial/LadrunoBondSlip.*",
    "SRC/material/nD/LogStrainNDMaterial.*", "SRC/material/nD/LogStrainKernel.h",
    "SRC/material/nD/FiniteStrainNDMaterial.h",
    "SRC/material/nD/LogStrain2D.*", "SRC/material/nD/FiniteStrainND2DMaterial.h",
    "SRC/material/nD/InitDefGradNDMaterial.*",
    "SRC/material/nD/StagedStrainNDMaterial.*",
    "SRC/analysis/integrator/CentralDifferenceLadruno.*",
    "SRC/analysis/integrator/CentralDifferenceSMS.*",
    "SRC/analysis/integrator/CentralDifferenceSMSConsistent.*",
    "SRC/analysis/integrator/LadrunoMassLumping.h",
    "SRC/analysis/integrator/LadrunoMassScaling.h",
    "SRC/analysis/integrator/LadrunoConsistentRefine.h",
    "SRC/analysis/integrator/LadrunoMassScalingEnergy.*",
    "SRC/analysis/integrator/LadrunoEnergyChannels.h",
    "SRC/analysis/integrator/LadrunoArcLength.*",
    "SRC/analysis/integrator/LadrunoDynamicRelaxation.*",
    "SRC/analysis/integrator/LadrunoFictitiousMass.h",
    "SRC/analysis/integrator/LadrunoIndirectControl.*",
    "SRC/analysis/integrator/ExplicitBathe.*",
    "SRC/recorder/LadrunoRecorder.*",
    "SRC/recorder/Ladruno_*.cpp", "SRC/recorder/Ladruno_*.h",
    "SRC/recorder/LadrunoMonitor*.cpp", "SRC/recorder/LadrunoMonitor*.h",
    "SRC/recorder/EnergyBalanceRecorder.*", "SRC/recorder/EnergyBalanceKernel.h",
    "SRC/convergenceTest/LadrunoStabilizedUnbalance.*",
    "SRC/analysis/handler/LadrunoProjectionHandler.*",
    "SRC/analysis/handler/LadrunoConstraintProjector.*",
    "SRC/analysis/handler/LadrunoProjectionConsumer.h",
    "SRC/analysis/handler/LadrunoContactHandler.*",
    "SRC/analysis/handler/LadrunoContactFE.*",
    "SRC/domain/constraints/LadrunoTie.*",
    "SRC/domain/contact/LadrunoContactDomain.*",
    "SRC/domain/contact/LadrunoContactSurface.*",
    "SRC/domain/contact/LadrunoContactKernel.h",
    "SRC/domain/contact/LadrunoContact2DKernel.h",   # ADR-85 T1a header; GLOBS debt discharged in T2
    "SRC/domain/contact/LadrunoContactProjection.h",
    "SRC/domain/contact/LadrunoFrictionKernel.h",
    "SRC/domain/contact/LadrunoMortarKernel.h",
    "SRC/domain/contact/LadrunoContactBucketSort.h",
    "SRC/domain/contact/LadrunoContactReemit.h",
    "SRC/domain/contact/LadrunoContactNormalField.h",
    "SRC/domain/contact/LadrunoEdgeKernel.h",
    "SRC/domain/pattern/ladrunoPorousOverlay/*.cpp", "SRC/domain/pattern/ladrunoPorousOverlay/*.h",
    "SRC/material/section/LadrunoShellModifierSection.*",
    # WP-122: stamped files that were missing here (a re-stamp would not have maintained them)
    "SRC/analysis/handler/LadrunoAutoPenaltyReduce.*",
    "SRC/analysis/handler/LadrunoContactAbort.*",
    "SRC/analysis/integrator/LadrunoSolverQuery.h",
    "SRC/material/nD/LadrunoCohesiveHingeBiaxial.*",
    "SRC/material/uniaxial/LadrunoCohesiveHinge.*",
    "SRC/utility/LadrunoParallelBuild.cpp",
    # WP-122: fork-added files that were never stamped, so the quirk lint never scanned
    # them (WP-120 R1: added on ladruno, absent upstream and at the merge-base)
    "SRC/Ladruno_mutation.h",
    "SRC/analysis/integrator/CriticalTimeStep.*",
    "SRC/analysis/integrator/LadrunoGeneralizedAlpha.*",
    "SRC/analysis/integrator/LadrunoHHT.*",
    "SRC/analysis/integrator/LadrunoLoadControl.*",
    "SRC/analysis/numberer/LadrunoParallelNumberer.*",
    "SRC/domain/pattern/drm/DRMHigherOrderNode.h",
    "SRC/element/LadrunoMassCache.h",
    "SRC/element/LadrunoResponseTokens.h",
    "SRC/interpreter/PythonMPIModule.cpp",
    "SRC/material/nD/LadrunoDamage.h",
    "SRC/utility/LadrunoThreads.*",
    "SRC/material/nD/ASDPlasticMaterial3D/HoekBrown_*.h",
    "SRC/material/nD/ASDPlasticMaterial3D/StiffSoil_*.h",
    "SRC/material/nD/ASDPlasticMaterial3D/test_HoekBrown.cpp",
    "SRC/material/nD/ASDPlasticMaterial3D/ElasticityModels/StiffSoil_EL.h",
    "SRC/material/nD/ASDPlasticMaterial3D/PlasticFlowDirections/HoekBrown_PF.h",
    "SRC/material/nD/ASDPlasticMaterial3D/PlasticFlowDirections/MohrCoulombTensionCutoff_PF.h",
    "SRC/material/nD/ASDPlasticMaterial3D/PlasticFlowDirections/StiffSoilCap_PF.h",
    "SRC/material/nD/ASDPlasticMaterial3D/PlasticFlowDirections/StiffSoilShear_PF.h",
    "SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/HoekBrown_YF.h",
    "SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/MohrCoulombTensionCutoff_YF.h",
    "SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/StiffSoilCap_YF.h",
    "SRC/material/nD/ASDPlasticMaterial3D/YieldFunctions/StiffSoilShear_YF.h",
]
SUFFIXES = {".cpp", ".h", ".hpp", ".cc", ".cxx"}


def authored_files() -> list[Path]:
    seen: dict[Path, None] = {}
    for g in GLOBS:
        for p in ROOT.glob(g):
            if p.is_file() and p.suffix in SUFFIXES:
                seen[p.resolve()] = None
    return sorted(seen)


# WP-122: --check only sees files already in GLOBS, so a fork file nobody added
# there stays unstamped -- and invisible to the quirk lint -- with a green check
# (how 31 files escaped, WP-120 R1). Any SRC source whose path names Ladruno is
# fork-authored by construction, so it must be in GLOBS. Not exhaustive: a fork
# file with a neutral name (CriticalTimeStep.cpp, the ASDPlastic kit headers)
# still relies on its author adding it to GLOBS.
LADRUNO_PATH = re.compile(r"ladruno", re.I)


def ladruno_named_outside_globs(root: Path, authored) -> list[Path]:
    auth = {p.resolve() for p in authored}
    return [p for p in sorted((root / "SRC").rglob("*"))
            if p.is_file() and p.suffix in SUFFIXES
            and LADRUNO_PATH.search(p.relative_to(root).as_posix())
            and p.resolve() not in auth]


def build_block(eol: str) -> str:
    art = (ROOT / "Ladruno_scripts" / "banner_ASCII.txt").read_text(
        encoding="utf-8").rstrip("\n").split("\n")
    lines = [START, RULE, "//"]
    for a in art:
        lines.append(("//  " + a).rstrip())
    lines.append("//")
    for c in CREDIT:
        lines.append("//  " + c)
    lines.append("//")
    lines.append("// Header auto-stamped by Ladruno_scripts/stamp_headers.py "
                 "(art: banner_ASCII.txt).")
    lines.append("// Do not hand-edit between the markers; edit the script/art "
                 "and re-run instead.")
    lines.append(RULE)
    lines.append(END)
    return eol.join(lines) + eol


_BLOCK_RE = re.compile(
    re.escape(START) + r".*?" + re.escape(END) + r"(?:\r?\n)?", re.S)


def restamp(text: str, block: str, eol: str) -> str:
    """Return text with the header block inserted or replaced (idempotent)."""
    if START in text:
        return _BLOCK_RE.sub(lambda _: block, text, count=1)
    # No existing block: insert after the leading /* ... */ comment if present.
    stripped = text.lstrip()
    if stripped.startswith("/*"):
        close = text.find("*/")
        nl = text.find("\n", close)
        if nl == -1:
            return text + eol + eol + block
        head, tail = text[:nl + 1], text[nl + 1:]
        return head + eol + block + tail
    return block + eol + text


def main() -> int:
    check = "--check" in sys.argv[1:]
    files = authored_files()
    changed: list[Path] = []
    for p in files:
        with open(p, "r", encoding="utf-8", newline="") as fh:
            raw = fh.read()
        bom = ""
        if raw.startswith("﻿"):
            bom, raw = "﻿", raw[1:]
        eol = "\r\n" if "\r\n" in raw else "\n"
        new = bom + restamp(raw, build_block(eol), eol)
        if new != bom + raw:
            changed.append(p)
            if not check:
                with open(p, "w", encoding="utf-8", newline="") as fh:
                    fh.write(new)

    rel = lambda p: p.relative_to(ROOT).as_posix()
    outside = ladruno_named_outside_globs(ROOT, files)
    if outside:
        print("Ladruno-named sources NOT in GLOBS ({}) -- add each to GLOBS, then stamp:"
              .format(len(outside)))
        for p in outside:
            print("  " + rel(p))
    if check:
        if changed:
            print("STALE / unstamped ({}):".format(len(changed)))
            for p in changed:
                print("  " + rel(p))
        if changed or outside:
            return 1
        print("All {} authored files carry a current header; every Ladruno-named source "
              "is in GLOBS.".format(len(files)))
        return 0

    print("Stamped {} of {} authored files (already-current files unchanged):"
          .format(len(changed), len(files)))
    for p in changed:
        print("  " + rel(p))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
