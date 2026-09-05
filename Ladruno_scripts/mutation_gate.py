#!/usr/bin/env python3
"""ADR-87 D2 — mutation gates: make "the tests are green" falsifiable.

A mutation build sabotages one physics family (see SRC/Ladruno_mutation.h).
This driver runs a family's suite against a normal binary and against a mutant
one, and reports what fraction of the tests NOTICED. A test that passes while
the physics is deleted is not testing the physics; the list of such survivors
is the actionable output of this tool, more than the score itself.

    # 1. normal build, record the baseline
    python Ladruno_scripts/mutation_gate.py run --family CONTINUUM \\
        --out artifacts/base.json

    # 2. rebuild with -DLADRUNO_MUTATE_FAMILY=CONTINUUM -DLADRUNO_MUTATE_MODE=ZERO
    python Ladruno_scripts/mutation_gate.py run --family CONTINUUM \\
        --expect CONTINUUM=ZERO --out artifacts/zero.json

    # 3. score
    python Ladruno_scripts/mutation_gate.py score --baseline artifacts/base.json \\
        --mutant artifacts/zero.json

WHY THE --expect CHECK IS NOT OPTIONAL
    The gate's claim is "the tests went red BECAUSE the physics was deleted".
    If a mutant build silently failed and the harness re-ran the previous
    binary, every test would pass and the tool would report the exact opposite
    of the truth -- "this suite cannot detect deleted physics" -- with no error
    anywhere. Every run therefore asks the binary what it is (`ladrunoMutation`)
    and refuses to proceed on a mismatch. Same fail-loud idiom as
    `ladrunoBuild` for stale-.pyd detection.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TESTS = REPO / "tests"

# --- family registry -------------------------------------------------------
# Each family names the suite that is supposed to be its evidence, and the
# floor its ZERO-mode score must clear. Thresholds are deliberately per-family
# and are RATCHETS: raise one when the suite improves, never lower it to make a
# red gate green -- that would be the same false comfort D2 exists to remove.
# INVARIANT: a family's `paths` must select exactly the suite whose elements
# actually carry that family's LADRUNO_MUTATE_* call sites. Selecting a wider
# net silently deflates the score -- tests of an UNMUTATED element cannot fail,
# so they look like "tests that ignore the physics" when they are nothing of
# the kind. (Caught during WP-2: the first CONTINUUM glob pulled in
# LadrunoBrick20 tests before that element had its gates.)
FAMILIES = {
    # 3D hex continuum lane: LadrunoBrick (std/bbar/URI/physical/SSP/EAS,
    # -geom finite + hypo) and LadrunoBrick20. Both are gated.
    "CONTINUUM": {
        "paths": ["test_ladrunoBrick*.py", "test_brick_inertia_skip.py"],
        # ZERO deletes the internal force outright: a suite that checks
        # equilibrium, reactions or displacements must overwhelmingly notice.
        "min_score": {"ZERO": 0.60, "SCALE": 0.40, "IDENT": 0.0},
    },
    # ADR 91 shell stiffness-modifier decorator. Gated at the SECTION accessors
    # (getStressResultant / getSectionTangent / getInitialTangent), not at an
    # element: the whole feature is a constitutive transform, so that is its
    # only physics surface. NOTE the expected survivors under ZERO are the three
    # pure PARSE-TIME refusal tests (R1/R3/R4) -- they never run an analysis, so
    # they cannot notice deleted physics, and that is the correct answer rather
    # than a gap. Gated.
    "SHELLMOD": {
        "paths": ["test_ladrunoShellModifier*.py"],
        "min_score": {"ZERO": 0.60, "SCALE": 0.30, "IDENT": 0.0},
    },
    # --- WP-3 lanes: gates NOT yet inserted in these families' elements. -----
    # Registering the paths here is deliberate (it documents the intended
    # evidence suite), but running a gate before the call sites exist would
    # score 0.0 and mean nothing. mutation_gate refuses to run a family whose
    # elements report no mutation, so the mistake is loud rather than silent.
    "UP": {"paths": ["test_ladruno_up*.py"], "min_score": {"ZERO": 0.50, "SCALE": 0.30, "IDENT": 0.0}},
    "CONTACT": {"paths": ["test_adr39_contact*.py", "test_ladruno_tie*.py"], "min_score": {"ZERO": 0.50, "SCALE": 0.30, "IDENT": 0.0}},
    "EXPLICIT": {"paths": ["test_*explicit*.py"], "min_score": {"ZERO": 0.50, "SCALE": 0.30, "IDENT": 0.0}},
    "SANISAND": {"paths": ["test_*manzari*.py", "test_*sanisand*.py"], "min_score": {"ZERO": 0.50, "SCALE": 0.30, "IDENT": 0.0}},
}

# Families whose C++ call sites are live. WP-3 moves the rest across as it
# inserts their gates; `run` refuses an ungated family so a 0.0 score can never
# be mistaken for "this suite is worthless".
GATED_FAMILIES = {"CONTINUUM", "SHELLMOD"}

PROBE = (
    "import json,sys\n"
    "try:\n"
    "    import opensees as ops\n"
    "except ImportError:\n"
    "    import openseespy.opensees as ops\n"
    "out={}\n"
    "try:\n"
    "    out['mutation']=ops.ladrunoMutation()\n"
    "except AttributeError:\n"
    "    out['mutation']='<verb-missing>'\n"
    "try:\n"
    "    out['build']=ops.ladrunoBuild()\n"
    "except Exception:\n"
    "    out['build']='<unknown>'\n"
    "print('LADRUNO_PROBE '+json.dumps(out))\n"
)


def probe_binary(python: str, cwd: Path) -> dict:
    """Ask the binary what it is. Never trust a result set without this.

    Runs in `cwd` -- the directory the extension module lives in. CI copies
    opensees.so into tests/ and the suite is run from there (`cd tests` is the
    Zone-A idiom), so probing from the repo root would raise ModuleNotFoundError
    even on a perfectly good build. Measured the hard way: run 33365505674.
    """
    res = subprocess.run([python, "-c", PROBE], capture_output=True, text=True, cwd=str(cwd))
    for line in res.stdout.splitlines():
        if line.startswith("LADRUNO_PROBE "):
            return json.loads(line[len("LADRUNO_PROBE "):])
    raise SystemExit(
        "FATAL: could not probe the OpenSees module.\n"
        f"stdout:\n{res.stdout}\nstderr:\n{res.stderr}"
    )


def resolve_paths(family: str) -> list[str]:
    """Test files as names RELATIVE to tests/ (pytest is run from there)."""
    spec = FAMILIES[family]
    found: list[str] = []
    for pat in spec["paths"]:
        found.extend(sorted(p.name for p in TESTS.glob(pat)))
    if not found:
        raise SystemExit(f"FATAL: family {family} selected no test files ({spec['paths']}).")
    return found


def cmd_run(args: argparse.Namespace) -> int:
    family = args.family.upper()
    if family not in FAMILIES:
        raise SystemExit(f"FATAL: unknown family {family}; known: {', '.join(FAMILIES)}")

    if family not in GATED_FAMILIES and args.expect:
        raise SystemExit(
            f"FATAL: family {family} has no LADRUNO_MUTATE_* call sites yet, so a\n"
            "mutation build cannot change its results and the score would be a\n"
            "meaningless 0.0. Insert the gates in that family's Element API\n"
            "accessors first (see SRC/Ladruno_mutation.h), then add it to\n"
            "GATED_FAMILIES."
        )

    python = args.python or sys.executable
    module_dir = Path(args.module_dir) if args.module_dir else TESTS
    info = probe_binary(python, module_dir)
    actual = info["mutation"]

    if actual == "<verb-missing>":
        raise SystemExit(
            "FATAL: this binary has no `ladrunoMutation` verb, so it predates the\n"
            "ADR-87 D2 gates and CANNOT be scored -- a green result here would be\n"
            "meaningless. Rebuild from a tree that contains SRC/Ladruno_mutation.h."
        )

    expected = args.expect or "none"
    if actual != expected:
        raise SystemExit(
            f"FATAL: binary identity mismatch. Expected mutation '{expected}', binary reports '{actual}'.\n"
            "Refusing to run: this is exactly the silent-stale-binary case that would\n"
            "make the gate report the OPPOSITE of the truth. Check that the build\n"
            "actually succeeded and that the interpreter is loading it (see the\n"
            "installed-.pth hijack note in Ladruno_internal/BUILD_GOTCHAS.md)."
        )

    paths = resolve_paths(family)
    with tempfile.TemporaryDirectory() as td:
        xml_path = Path(td) / "results.xml"
        cmd = [python, "-m", "pytest", "-p", "no:randomly", "--tb=no", "-q",
               f"--junit-xml={xml_path}", *paths]
        if args.marker:
            cmd.extend(["-m", args.marker])
        print("running:", " ".join(cmd[:6]), f"... ({len(paths)} files)", flush=True)
        subprocess.run(cmd, cwd=str(module_dir))   # exit code ignored: red is expected under mutation

        if not xml_path.exists():
            raise SystemExit("FATAL: pytest produced no JUnit XML -- the run did not start.")
        results = parse_junit(xml_path)

    payload = {
        "family": family,
        "mutation": actual,
        "build": info.get("build", "<unknown>"),
        "n_tests": len(results),
        "results": results,
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    n_pass = sum(1 for v in results.values() if v == "passed")
    print(f"{family}: {n_pass}/{len(results)} passed  (mutation={actual})  -> {args.out}")
    return 0


def parse_junit(path: Path) -> dict:
    """testcase -> passed | failed | skipped."""
    root = ET.parse(path).getroot()
    out: dict[str, str] = {}
    for tc in root.iter("testcase"):
        name = f"{tc.get('classname', '')}::{tc.get('name', '')}"
        status = "passed"
        for child in tc:
            tag = child.tag.lower()
            if tag in ("failure", "error"):
                status = "failed"
                break
            if tag == "skipped":
                status = "skipped"
                break
        out[name] = status
    return out


def cmd_score(args: argparse.Namespace) -> int:
    base = json.loads(Path(args.baseline).read_text(encoding="utf-8"))
    mut = json.loads(Path(args.mutant).read_text(encoding="utf-8"))

    if base["mutation"] != "none":
        raise SystemExit(f"FATAL: baseline was itself a mutant build ({base['mutation']}).")
    if mut["mutation"] == "none":
        raise SystemExit(
            "FATAL: the 'mutant' run reports mutation='none' -- nothing was sabotaged,\n"
            "so a green suite proves nothing. The mutation build did not take effect."
        )
    if base["family"] != mut["family"]:
        raise SystemExit(f"FATAL: family mismatch ({base['family']} vs {mut['family']}).")

    mode = mut["mutation"].split("=")[-1]

    # Only tests that PASSED in the baseline can detect anything. A test already
    # failing (or skipped for a missing optional dep) carries no information, so
    # including it would silently deflate or inflate the score.
    detectors = [t for t, st in base["results"].items() if st == "passed"]
    if not detectors:
        raise SystemExit("FATAL: no baseline-passing tests -- nothing could detect a mutation.")

    killed, survived, vanished = [], [], []
    for t in detectors:
        st = mut["results"].get(t)
        if st is None:
            vanished.append(t)          # collection differed between runs
        elif st == "failed":
            killed.append(t)
        else:
            survived.append(t)

    score = len(killed) / len(detectors)
    floor = args.min_score
    if floor is None:
        floor = FAMILIES.get(mut["family"], {}).get("min_score", {}).get(mode, 0.0)

    print(f"\nADR-87 D2 mutation gate - {mut['family']} / {mode}")
    print(f"  baseline passing tests : {len(detectors)}")
    print(f"  killed by the mutation : {len(killed)}")
    print(f"  SURVIVED (do not test the physics): {len(survived)}")
    if vanished:
        print(f"  not collected in mutant run: {len(vanished)}")
    print(f"  mutation score : {score:.3f}   (floor {floor:.3f})")

    if survived:
        print("\n  survivors - each is a test that passes with the physics deleted:")
        for t in sorted(survived)[: args.list_limit]:
            print(f"    {t}")
        if len(survived) > args.list_limit:
            print(f"    ... and {len(survived) - args.list_limit} more "
                  f"(raise --list-limit to see them)")

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps({
            "family": mut["family"], "mode": mode, "score": score, "floor": floor,
            "n_detectors": len(detectors), "killed": sorted(killed),
            "survived": sorted(survived), "not_collected": sorted(vanished),
        }, indent=2), encoding="utf-8")

    if score < floor:
        print(f"\nGATE FAILED: score {score:.3f} < floor {floor:.3f}.")
        print("The suite does not adequately detect deleted physics. Fix the TESTS,")
        print("not the floor -- lowering the floor to pass is the failure mode this")
        print("gate exists to prevent.")
        return 1
    print("\nGATE PASSED.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="ADR-87 D2 mutation gate driver")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run a family's suite and record per-test results")
    r.add_argument("--family", required=True)
    r.add_argument("--expect", help="mutation string the binary must report (default 'none')")
    r.add_argument("--out", required=True)
    r.add_argument("--python", help="interpreter to use (default: this one)")
    r.add_argument("--marker", help="extra pytest -m expression, e.g. zone_a")
    r.add_argument("--module-dir",
                   help="directory holding the built opensees module and the test files "
                        "(default: tests/). Both the probe and pytest run here.")
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("score", help="compare a baseline and a mutant run")
    s.add_argument("--baseline", required=True)
    s.add_argument("--mutant", required=True)
    s.add_argument("--min-score", type=float, default=None)
    s.add_argument("--out")
    s.add_argument("--list-limit", type=int, default=25)
    s.set_defaults(func=cmd_score)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
