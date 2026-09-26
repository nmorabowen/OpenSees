#!/usr/bin/env python3
"""WP-123 mutation rows: prove tests/test_ladruno_undamped_couplings.py fails on one-line
breaks of the shared code, and record which breaks are EQUIVALENT (undetectable) and why.

Each row edits the sources (exact, count-checked replacements), rebuilds ONLY opensees.pyd
(`build.bat OpenSeesPy`), runs the new test and the four element batteries as SEPARATE
processes (a crash in one must not hide the others), then restores every edited file in a
`finally`. After the rows, run a full 5-target build.bat (a named target leaves the other
dist binaries stale — AGENTS.md rule 2).

    <py3.12> mutation_rows.py [ROW ...]        (default: B C D E)

Needs cmd.exe + the oneAPI env (Ladruno_scripts/setup_env.bat). Writes a summary table.
"""
from __future__ import annotations

import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PY = r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe"
HDR = ROOT / "SRC/element/ladrunoEmbeddedRebar/LadrunoUndampedElement.h"
ELEM_CPP = ROOT / "SRC/element/Element.cpp"

GETDAMP_BLOCK = re.compile(
    r"  const Matrix &getDamp\(void\)\r?\n  \{\r?\n.*?return zeroDamp;\r?\n  \}\r?\n", re.S)
SETR_NOOP = "int setRayleighDampingFactors(double, double, double, double) { return 0; }"
ZERO_RET = re.compile(r"(    zeroDamp\.Zero\(\);\r?\n)(    return zeroDamp;)")
# the getDamp self-heal in vanilla Element.cpp (the 2026-07-28 fork fix), and ONLY that site
ELEM_GETDAMP_HEAL = re.compile(
    r"(Element::getDamp\(void\)\s*\r?\n\{\r?\n\s*if \(index\s+== -1\) \{\r?\n\s*)"
    r"this->Element::setRayleighDampingFactors\(alphaM, betaK, betaK0, betaKc\);(\s*// Ladruno)")

ROWS = {
    "B": ("base STORES the Rayleigh factors (setRayleighDampingFactors forwards to Element)",
          [(HDR, "lit", SETR_NOOP,
            "int setRayleighDampingFactors(double a, double b, double c, double d) "
            "{ return this->Element::setRayleighDampingFactors(a, b, c, d); }   // MUTATION B")]),
    "C": ("getDamp returns a NONZERO entry (C(0,0) = 1e3)",
          [(HDR, "re", ZERO_RET, r"\1    if (n > 0) zeroDamp(0, 0) = 1.0e3;   // MUTATION C\n\2")]),
    "D": ("no getDamp override AND Element.cpp getDamp self-heal un-qualified (the #219 shape)",
          [(HDR, "re", GETDAMP_BLOCK, "  // MUTATION D: getDamp override removed\n"),
           (ELEM_CPP, "re", ELEM_GETDAMP_HEAL,
            r"\1this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);\2 MUTATION D")]),
    "E": ("no getDamp override (base getDamp, Element.cpp fix intact) -- expected EQUIVALENT",
          [(HDR, "re", GETDAMP_BLOCK, "  // MUTATION E: getDamp override removed\n")]),
}

TESTS = [
    "tests/test_ladruno_undamped_couplings.py",
    "tests/test_ladrunoDistributingCoupling_element.py",
    "tests/test_ladrunoKinematicCoupling_element.py",
    "tests/test_ladrunoEmbeddedNode_element.py",
    "tests/test_ladrunoEmbeddedRebar_element.py",
]


def apply(edits):
    saved = {}
    for path, kind, old, new in edits:
        raw = path.read_bytes()
        saved.setdefault(path, raw)
        text = raw.decode("utf-8")
        if kind == "lit":
            assert text.count(old) == 1, (path, old)
            text = text.replace(old, new)
        else:
            text, n = old.subn(new, text)
            assert n == 1, (path, old.pattern, n)
        path.write_bytes(text.encode("utf-8"))
    return saved


def build_pyd(log):
    t0 = time.time()
    r = subprocess.run('cmd /c "call Ladruno_scripts\\setup_env.bat && Ladruno_scripts\\build.bat OpenSeesPy"',
                       cwd=ROOT, shell=True, capture_output=True, text=True, errors="replace")
    log.write_text(r.stdout + r.stderr, encoding="utf-8")
    return r.returncode, time.time() - t0


def run_tests():
    out = []
    for t in TESTS:
        r = subprocess.run([PY, "-S", str(ROOT / "Ladruno_implementation/wp123_undamped/run_pytest.py"), t],
                           cwd=ROOT, capture_output=True, text=True, errors="replace")
        tail = [ln for ln in r.stdout.strip().splitlines() if ln.strip()]
        summary = tail[-1] if tail else "(no output)"
        if r.returncode not in (0, 1):
            summary = f"CRASHED rc={r.returncode} (0x{r.returncode & 0xFFFFFFFF:08X}) | {summary}"
        out.append((Path(t).name, r.returncode, summary))
    return out


def main():
    rows = sys.argv[1:] or ["B", "C", "D", "E"]
    logdir = ROOT / "build" / "wp123_mutation_logs"
    logdir.mkdir(parents=True, exist_ok=True)
    report = []
    for row in rows:
        desc, edits = ROWS[row]
        print(f"=== row {row}: {desc}", flush=True)
        saved = apply(edits)
        try:
            rc, secs = build_pyd(logdir / f"build_{row}.log")
            if rc != 0:
                report.append((row, desc, [("BUILD FAILED", rc, f"see {logdir / f'build_{row}.log'}")]))
                print(f"  build FAILED rc={rc}", flush=True)
                continue
            print(f"  built opensees.pyd in {secs:.0f} s", flush=True)
            res = run_tests()
            for name, rc2, summ in res:
                print(f"  {name:48s} rc={rc2}  {summ}", flush=True)
            report.append((row, desc, res))
        finally:
            for path, raw in saved.items():
                path.write_bytes(raw)
            print("  sources restored", flush=True)
    print("\n| Row | Mutation | " + " | ".join(Path(t).stem.replace("test_", "") for t in TESTS) + " |")
    print("|---|---|" + "---|" * len(TESTS))
    for row, desc, res in report:
        print(f"| {row} | {desc} | " + " | ".join(s for _, _, s in res) + " |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
