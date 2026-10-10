#!/usr/bin/env python3
"""MKL-compatibility gate (WP-148). The MKL-gated sources must build against the
OLDEST MKL a fork build uses: esmeralda's spack oneMKL 2024.2.2.

Why: WP-132 (#864) used `MKL_CBWR_AVX10`, a oneMKL 2025.0 macro. The Windows box
has 2025.1 and built; esmeralda has 2024.2.2, so the Linux PARDISO opt-in
(`-DLADRUNO_MKL_PARDISO_LINUX=ON`) stopped compiling and nobody noticed for two
days (fixed in #886). No other job compiles these sources on Linux: Zone-A
configures without MKL, so `SRC/system_of_eqn/linearSOE/pardiso/`, the FEAST code
under `_LADRUNO_MKL_FEAST` and the `_PARDISO` blocks never reach a compiler there.
LEDGER_quirks: "An MKL symbol newer than esmeralda's oneMKL 2024.2 compiles on
Windows and breaks the Linux PARDISO build".

Input: a build directory configured with CMAKE_EXPORT_COMPILE_COMMANDS=ON and both
Linux opt-ins (LADRUNO_MKL_PARDISO_LINUX, LADRUNO_MKL_FEAST_LINUX) against the
oneMKL 2024.2.2 headers + libraries. Nothing is built beforehand.

  M1 compile   Every translation unit that touches MKL compiles, with the exact
               command CMake generated for it. "Touches MKL" = includes an MKL
               header (directly, or through ProfilerRunMeta.h, which includes
               mkl_service.h under _PARDISO) or declares an MKL routine itself.
               Catches a too-new MKL macro, type or header-declared function.
  M2 symbols   Every MKL symbol those objects reference is defined by the MKL
               libraries the opt-in links (LADRUNO_MKL_LP64/SEQ/CORE from
               CMakeCache.txt). Needed because the FEAST sources declare
               dfeast_srci / feastinit / pardiso in their own extern "C" blocks,
               so a too-new routine there compiles and only fails at link time.

The gate refuses to pass on nothing (exit 2): the configured MKL headers must be
the expected version, and each REQUIRED unit must be compiled with the macro that
switches its MKL code on. `--self-test` also proves both rules still fire, on two
canaries: M1 on the incident itself (`MKL_CBWR_AVX10` used without a guard), M2 on
a declared routine no MKL defines.

Exit codes: 0 clean, 1 finding, 2 the gate cannot evaluate (no compile database,
a required unit missing or without its macro, wrong MKL version, a self-test
canary that did not fire).

    python ci/check_mkl_compat.py --build-dir build/mkl [--self-test] [-j N]
"""
import argparse
import os
import re
import shlex
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

OLDEST_MKL = 20240002  # esmeralda: spack intel-oneapi-mkl 2024.2.2

INCLUDES_MKL = re.compile(
    r'#\s*include\s*[<"](?:[\w./]*/)?(?:mkl\w*\.h|ProfilerRunMeta\.h)[>"]')
DECLARES_MKL = re.compile(
    r'\b(?:void|int|double|MKL_INT)\s+\*?\s*'
    r'(?:pardiso\w*|[sdcz]?feast\w*|mkl_\w+|MKL_\w+)\s*\(')
# C-API names MKL must provide. Service routines are exported under their MKL_
# spelling (mkl_service.h: `#define mkl_get_max_threads MKL_Get_Max_Threads`).
MKL_SYMBOL = re.compile(
    r'^(?:mkl_|MKL_|pardiso|PARDISO|[sdcz]feast_|feastinit|'
    r'cluster_sparse_solver|LAPACKE_|cblas_)')

# Units whose MKL path must really be compiled, and the macro that switches it
# on (None: the unit is compiled only when MKL is on, so it has no switch).
REQUIRED = {
    "SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSolver.cpp": None,
    "SRC/system_of_eqn/eigenSOE/FeastEigenSolver.cpp": "_LADRUNO_MKL_FEAST",
    "SRC/system_of_eqn/eigenSOE/LadrunoBlockZKernel.cpp": "_LADRUNO_MKL_FEAST",
    "SRC/tcl/commands.cpp": "_PARDISO",
    "SRC/interpreter/OpenSeesCommands.cpp": "_PARDISO",
}
MKL_LIB_VARS = ("LADRUNO_MKL_LP64", "LADRUNO_MKL_SEQ", "LADRUNO_MKL_CORE")

CANARY_M1 = """\
// WP-148 self-test, M1: the #864 incident. MKL_CBWR_AVX10 is oneMKL 2025.0+,
// so against 2024.2 headers this must NOT compile.
#include <mkl_types.h>
int ladruno_mkl_canary_m1 = MKL_CBWR_AVX10;
"""
CANARY_M2_SYMBOL = "mkl_ladruno_canary_absent"
CANARY_M2 = """\
// WP-148 self-test, M2: a routine declared by hand (the FEAST pattern) that no
// MKL defines. It compiles; the symbol rule must report it.
extern "C" void %s(void);
void ladruno_mkl_canary_m2(void) { %s(); }
""" % (CANARY_M2_SYMBOL, CANARY_M2_SYMBOL)


class CannotEvaluate(Exception):
    pass


def read_cache(build):
    cache = {}
    path = build / "CMakeCache.txt"
    if not path.is_file():
        raise CannotEvaluate(f"{path} not found -- configure the build directory first")
    for line in path.read_text(errors="replace").splitlines():
        m = re.match(r"^([A-Za-z_0-9]+):[A-Z_]+=(.*)$", line)
        if m:
            cache[m.group(1)] = m.group(2)
    return cache


def mkl_version(include_dir):
    hdr = Path(include_dir) / "mkl_version.h"
    if not hdr.is_file():
        raise CannotEvaluate(f"{hdr} not found (LADRUNO_MKL_INCLUDE={include_dir})")
    m = re.search(r"#define\s+INTEL_MKL_VERSION\s+(\d+)", hdr.read_text(errors="replace"))
    if not m:
        raise CannotEvaluate(f"no INTEL_MKL_VERSION in {hdr}")
    return int(m.group(1))


def load_units(build, root):
    import json
    db = build / "compile_commands.json"
    if not db.is_file():
        raise CannotEvaluate(f"{db} not found -- configure with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON")
    units = []
    for e in json.loads(db.read_text()):
        src = Path(e["file"])
        if not src.is_absolute():
            src = Path(e["directory"]) / src
        try:
            rel = src.resolve().relative_to(root).as_posix()
        except ValueError:
            continue
        if not rel.startswith("SRC/"):
            continue  # OTHER/ third-party code is not ours to gate
        text = src.read_text(errors="replace")
        if not (INCLUDES_MKL.search(text) or DECLARES_MKL.search(text)):
            continue
        args = e["arguments"] if "arguments" in e else shlex.split(e["command"])
        units.append({"rel": rel, "dir": e["directory"], "args": args})
    return units


def output_of(args):
    for i, a in enumerate(args):
        if a == "-o" and i + 1 < len(args):
            return args[i + 1]
    raise CannotEvaluate("compile command has no -o: " + " ".join(args[:3]) + " ...")


def has_macro(args, name):
    return any(a == "-D" + name or a.startswith("-D" + name + "=") for a in args)


def check_required(units):
    problems = []
    for rel, macro in REQUIRED.items():
        mine = [u for u in units if u["rel"] == rel]
        if not mine:
            problems.append(f"{rel} is not in the compile database (or no longer "
                            "touches MKL) -- the gate would not compile it")
        elif macro and not any(has_macro(u["args"], macro) for u in mine):
            problems.append(f"{rel} is compiled without -D{macro}, so its MKL code "
                            "is preprocessed away -- is the opt-in on?")
    if problems:
        raise CannotEvaluate("\n  ".join(["required units:"] + problems))


def compile_unit(u):
    out = Path(u["dir"]) / output_of(u["args"])
    out.parent.mkdir(parents=True, exist_ok=True)
    p = subprocess.run(u["args"], cwd=u["dir"], capture_output=True, text=True)
    return u, out, p.returncode, p.stdout + p.stderr


def nm_undefined(obj):
    p = subprocess.run(["nm", "-u", "--format=posix", str(obj)],
                       capture_output=True, text=True, check=True)
    return {ln.split()[0].split("@")[0] for ln in p.stdout.splitlines() if ln.strip()}


def nm_defined(libs):
    names = set()
    for lib in libs:
        p = subprocess.run(["nm", "-D", "--defined-only", "--format=posix", str(lib)],
                           capture_output=True, text=True, check=True)
        names |= {ln.split()[0].split("@")[0] for ln in p.stdout.splitlines() if ln.strip()}
    return names


def mkl_references(objects):
    """{MKL symbol: [objects that reference it]}"""
    refs = {}
    for obj in objects:
        for s in nm_undefined(obj):
            if MKL_SYMBOL.match(s):
                refs.setdefault(s, []).append(obj)
    return refs


def missing_mkl_symbols(objects, defined):
    return {s: o for s, o in mkl_references(objects).items() if s not in defined}


def canary_unit(template, text, tmp, name):
    """A copy of `template`'s compile command that compiles `text` instead."""
    src = Path(tmp) / f"{name}.cpp"
    src.write_text(text)
    # CMake passes the source as an absolute path into the repo. The dependency
    # flags go: they would overwrite the real unit's .d file in the build tree.
    args, skip = [], False
    for a in template["args"]:
        if skip:
            skip = False
        elif a in ("-MF", "-MT", "-MQ"):
            skip = True
        elif a not in ("-MD", "-MMD"):
            args.append(str(src) if a.endswith("/" + template["rel"]) else a)
    i = args.index("-o")
    args[i + 1] = str(Path(tmp) / f"{name}.o")
    if str(src) not in args:
        raise CannotEvaluate(f"cannot locate the source argument in {template['rel']}'s command")
    return {"rel": f"<{name}>", "dir": template["dir"], "args": args}


def self_test(template, defined, tmp):
    failures = []
    _, _, rc, log = compile_unit(canary_unit(template, CANARY_M1, tmp, "canary_m1"))
    if rc == 0:
        failures.append("M1 canary (MKL_CBWR_AVX10 unguarded) COMPILED -- these headers "
                        "are newer than 2024.2 or not the ones in use; M1 has no teeth")
    elif "MKL_CBWR_AVX10" not in log:
        failures.append("M1 canary failed for a reason other than MKL_CBWR_AVX10:\n" + log[-1500:])
    _, obj, rc, log = compile_unit(canary_unit(template, CANARY_M2, tmp, "canary_m2"))
    if rc != 0:
        failures.append("M2 canary did not compile:\n" + log[-1500:])
    elif CANARY_M2_SYMBOL not in missing_mkl_symbols([obj], defined):
        failures.append(f"M2 canary: {CANARY_M2_SYMBOL} was NOT reported missing -- M2 has no teeth")
    return failures


def main():
    ap = argparse.ArgumentParser(description="MKL-compatibility gate (WP-148).")
    ap.add_argument("--build-dir", type=Path, required=True)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--expect-mkl-version", type=int, default=OLDEST_MKL)
    ap.add_argument("--self-test", action="store_true",
                    help="also prove M1 and M2 fire, on two canaries")
    ap.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 2)
    a = ap.parse_args()
    root, build = a.root.resolve(), a.build_dir.resolve()
    try:
        cache = read_cache(build)
        inc = cache.get("LADRUNO_MKL_INCLUDE", "")
        if not inc or inc.endswith("-NOTFOUND"):
            raise CannotEvaluate("LADRUNO_MKL_INCLUDE is not set -- configure with "
                                 "-DLADRUNO_MKL_PARDISO_LINUX=ON -DMKL_RT_HINT=<mkl>/lib")
        ver = mkl_version(inc)
        if ver != a.expect_mkl_version:
            raise CannotEvaluate(f"configured MKL is INTEL_MKL_VERSION {ver}, expected "
                                 f"{a.expect_mkl_version} (the oldest MKL a fork build uses)")
        libs = [cache.get(v, "") for v in MKL_LIB_VARS]
        if not all(libs) or any(l.endswith("-NOTFOUND") for l in libs):
            raise CannotEvaluate(f"{'/'.join(MKL_LIB_VARS)} not all set in CMakeCache.txt")
        units = load_units(build, root)
        check_required(units)
    except CannotEvaluate as e:
        print(f"check_mkl_compat: cannot evaluate -- {e}")
        return 2

    print(f"check_mkl_compat: MKL {ver} ({inc}); {len(units)} unit(s) touch MKL:")
    for u in units:
        print(f"  {u['rel']}")
    findings = []
    with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as pool:
        results = list(pool.map(compile_unit, units))
    objects = []
    for u, out, rc, log in results:
        if rc != 0:
            errs = [ln for ln in log.splitlines() if "error" in ln]
            findings.append(f"M1 {u['rel']}: does not compile against MKL {ver}:\n    "
                            + "\n    ".join(errs[:8] or log.splitlines()[-8:]))
        else:
            objects.append(out)
    defined = nm_defined(libs)
    referenced = mkl_references(objects)
    for sym, objs in sorted(referenced.items()):
        if sym in defined:
            continue
        names = sorted({Path(o).name for o in objs})
        findings.append(f"M2 {sym}: referenced by {', '.join(names)} but not defined by "
                        f"MKL {ver} ({', '.join(Path(l).name for l in libs)})")

    if a.self_test:
        template = next(u for u in units if u["rel"] == next(iter(REQUIRED)))
        with tempfile.TemporaryDirectory() as tmp:
            failures = self_test(template, defined, tmp)
        if failures:
            for f in failures:
                print("SELF-TEST " + f)
            print("check_mkl_compat: self-test FAILED -- the gate cannot be trusted")
            return 2
        print("check_mkl_compat: self-test OK (M1 and M2 canaries both fire)")

    for f in findings:
        print(f)
    print(f"check_mkl_compat: MKL symbols referenced: {' '.join(sorted(referenced)) or '-'}")
    print(f"check_mkl_compat: {len(objects)}/{len(units)} unit(s) compiled, "
          f"{len(referenced)} MKL symbol(s) referenced, {len(findings)} finding(s)")
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
