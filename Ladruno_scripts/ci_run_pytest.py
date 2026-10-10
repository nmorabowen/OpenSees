"""Run pytest against ONE known engine build, and prove it (WP-179).

    py -3.12 -S Ladruno_scripts/ci_run_pytest.py [--dist DIR] [--expect-sha SHA] -- <pytest args>

The Windows-only gates in one command (BUILD_GOTCHAS section 4c):

    py -3.12 -S Ladruno_scripts/ci_run_pytest.py -- $(py -3.12 ci/check_quirk_patterns.py --list-platform-tests)

Every way local Windows test evidence went wrong in WP-175..178 is closed here:

* `-S` is REQUIRED. Without it `site` runs a boot `.pth` that imports
  `opensees` from whatever worktree it was last wired to, before `sys.path` is
  ours (WP-175/176). `site-packages` is then added by path, so pytest and numpy
  resolve, but no `.pth` IMPORT line runs.
* `dist\\bin` is pinned for this process (`sys.path` + `add_dll_directory`)
  and for children (`PYTHONPATH`), and `opensees.__file__` is asserted.
* `--expect-sha` (default `$GITHUB_SHA` when set, else no check) must equal
  `opensees.ladrunoBuild()`, which proves the tests run on the commit just built.
  The stamp is taken at CMake CONFIGURE time (BUILD_GOTCHAS "ladrunoBuild lags"),
  so touch CMakeLists.txt before the build when you pass it.

Exit codes: pytest's own (0-5), or 90 = wrong/stale engine, 91 = not started with -S.
Deliberately outside pytest's range, so a caller can tell "the tests ran and
failed" from "the tests never ran".
"""
import argparse
import os
import sys
import sysconfig

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def main(argv):
    if not sys.flags.no_site:
        print("ci_run_pytest: run me with `python -S` -- without it a site-packages boot .pth "
              "can load another build's opensees before sys.path is pinned (WP-175/176).",
              file=sys.stderr)
        return 91
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("--dist", default=os.path.join(ROOT, "dist", "bin"))
    ap.add_argument("--expect-sha", default=os.environ.get("GITHUB_SHA", ""))
    ours, pytest_args = (argv[:argv.index("--")], argv[argv.index("--") + 1:]) if "--" in argv \
        else (argv, [])
    args = ap.parse_args(ours)

    dist = os.path.abspath(args.dist)
    tests = os.path.join(ROOT, "tests")
    site_pkgs = sysconfig.get_paths()["purelib"]
    os.add_dll_directory(dist)
    sys.path[:0] = [dist, tests]
    sys.path.append(site_pkgs)
    os.environ["PYTHONPATH"] = os.pathsep.join([dist, tests, site_pkgs])
    os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

    if "opensees" in sys.modules:
        print("ci_run_pytest: opensees was imported before the paths were pinned", file=sys.stderr)
        return 90
    import opensees
    want = os.path.normcase(os.path.join(dist, "opensees.pyd" if os.name == "nt" else "opensees.so"))
    got = os.path.normcase(os.path.abspath(opensees.__file__))
    build = str(opensees.ladrunoBuild())
    print("ci_run_pytest: engine %s  ladrunoBuild %s" % (opensees.__file__, build), flush=True)
    if got != want:
        print("ci_run_pytest: loaded %s, expected %s" % (got, want), file=sys.stderr)
        return 90
    if args.expect_sha and not build.startswith(args.expect_sha) and not args.expect_sha.startswith(build):
        print("ci_run_pytest: engine was built from %s, expected %s -- dist\\bin is stale or the "
              "stamp lags (force a CMake reconfigure before building)" % (build, args.expect_sha),
              file=sys.stderr)
        return 90

    os.chdir(tests)
    import pytest
    return pytest.main(pytest_args + ["-p", "no:cacheprovider"])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
