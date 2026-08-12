"""Driver for the `constraints Auto` cross-rank reduction test.

Runs the three serial references and the 2-rank deck, scrapes every `PVAL =`
line the handler's `-verbose` report prints, and asserts the four expectations.

    python run.py <build-dir>

<build-dir> is the directory holding OpenSees.exe / OpenSeesMP.exe -- normally
`build\\build\\Release` under the tree you built. The binaries are named
explicitly on purpose: a bare `OpenSees` on PATH has resolved to another
session's build on this machine more than once.
"""

import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PVAL = re.compile(r"PVAL\s*=\s*([0-9.eE+-]+)")

# part -> the PVAL that part's element population implies. See common.tcl.
SOFT, STIFF, WHOLE = 1.0e5, 1.0e9, 1.0e7


def oneapi_dirs():
    """The oneAPI directories OpenSeesMP.exe needs on PATH to even load.

    Without them every rank dies at load with exit status 0xC0000135
    (STATUS_DLL_NOT_FOUND) before printing anything -- which reads as "no PVAL
    lines" and is easy to mistake for a failed assertion rather than a broken
    launch. Resolving mpiexec.exe by absolute path is NOT enough; impi.dll is
    the problem.

    Built explicitly rather than by shelling through setup_env.bat, because
    that script inherits PATH and misbehaves when the parent shell is not
    cmd.exe (from Git Bash it fails to find its own tools and silently leaves
    mpiexec unresolvable, still exiting 0).
    """
    root = os.environ.get("I_MPI_ROOT")
    oneapi = os.environ.get("ONEAPI_ROOT", r"C:\Program Files (x86)\Intel\oneAPI")
    mpiroot = root or os.path.join(oneapi, "mpi", "latest")
    cands = [
        os.path.join(mpiroot, "bin"),
        # libfabric is what PMPI_Init actually needs, and it is NOT under
        # mpi\latest\libfabric -- Intel MPI's own vars.bat puts
        # %I_MPI_ROOT%\opt\mpi\libfabric\bin on PATH. Miss it and every rank
        # gets past DLL load only to die in MPIDI_OFI_mpi_init_hook with
        # "Unknown error class", which names nothing.
        os.path.join(mpiroot, "opt", "mpi", "libfabric", "bin"),
        os.path.join(oneapi, "compiler", "latest", "bin"),
    ]
    return [c for c in cands if c and os.path.isdir(c)]


def find_mpiexec(dirs):
    for d in dirs:
        c = os.path.join(d, "mpiexec.exe")
        if os.path.exists(c):
            return c
    return "mpiexec"


def run(cmd, env):
    p = subprocess.run(cmd, cwd=HERE, env=env, capture_output=True, text=True,
                       timeout=300)
    return p.stdout + p.stderr


def pvals(out):
    return [float(v) for v in PVAL.findall(out)]


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    bd = os.path.abspath(sys.argv[1])
    ops = os.path.join(bd, "OpenSees.exe")
    opsmp = os.path.join(bd, "OpenSeesMP.exe")
    for exe in (ops, opsmp):
        if not os.path.exists(exe):
            print("MISSING %s" % exe)
            return 2

    env = dict(os.environ)
    # find dist/lib/tcl8.6 by walking up -- bd may be build\build\Release or dist\bin
    tcllib, d = None, bd
    for _ in range(6):
        cand = os.path.join(d, "dist", "lib", "tcl8.6")
        if os.path.isdir(cand):
            tcllib = cand
            break
        nd = os.path.dirname(d)
        if nd == d:
            break
        d = nd
    if tcllib:
        env["TCL_LIBRARY"] = tcllib

    dirs = oneapi_dirs()
    # bd first: the MKL / libiomp DLLs sit next to the binaries in dist\bin
    env["PATH"] = os.pathsep.join([bd] + dirs + [env.get("PATH", "")])

    failures = []

    def check(name, got, want, n=1):
        ok = len(got) == n and all(abs(g - want) <= 1e-6 * want for g in got)
        print("%-28s %-10s got=%s want=%s x%d" %
              (name, "PASS" if ok else "FAIL",
               [("%.6e" % g) for g in got], "%.6e" % want, n))
        if not ok:
            failures.append(name)

    for part, want in (("soft", SOFT), ("stiff", STIFF), ("whole", WHOLE)):
        out = run([ops, "serial.tcl", part], env)
        check("serial/" + part, pvals(out), want)

    out = run([find_mpiexec(dirs), "-n", "2", opsmp, "mp.tcl"], env)
    got = pvals(out)
    if not got and "C0000135" in out.upper():
        print("  ^ every rank died at load (STATUS_DLL_NOT_FOUND), so this is a "
              "broken launch, not a failed assertion.")
    # Both ranks must report the whole-model value. Pre-fix this is
    # [1e5, 1e9] in some order.
    check("mp 2-rank (both ranks)", got, WHOLE, n=2)
    if sorted(got) == [SOFT, STIFF]:
        print("  ^ this is exactly the pre-fix signature: each rank sized its "
              "penalty from its own elements.")

    print()
    if failures:
        print("FAILED: " + ", ".join(failures))
        return 1
    print("ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
